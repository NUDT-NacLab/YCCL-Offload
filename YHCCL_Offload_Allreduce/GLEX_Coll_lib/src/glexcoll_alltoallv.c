#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <time.h>

#include "glexcoll.h"
#include "glexalltoall.h"

//整型数据压缩库
#include "vint.h"
#include "vsimple.h"
#include "vp4.h"
#include "bitpack.h"
#include "eliasfano.h"

static int group_rank;
static int group_procn;
MPI_Comm Comm_group;
static int internode_procn;
static int internode_rank;
static int flag_control = 0;
double *statistics_time;
static int statistics_count = 0;
int AlltoallvAlgorithm;
enum AlltoallvReorderTYPE
{
    reverse_order,
    random_order,
    diavation_order,
    ring_distance_order
};
int AlltoallvReorderType;
//alltoallv算法中同时允许发送的消息请求的数量
int alltoallV_ongoingmsgN = 2;
int *random_order_vec;
char *diavation_order_vec;
#define DIAVATION_VEC_ELEM_SIZE 8
void alltoallv_update_communication_order();
double *total_time_MAX, *total_time_MIN;

static MPI_Request *reqSV;
static MPI_Request *reqRV;
static MPI_Status *statusV;

//注册用于发送接收消息的缓冲区
char *AlltoallvBufSends[64];
char *AlltoallvBufRecvs[64];

//用于rdma远程内存访问
uint64_t *rmt_disps;

void GLEXCOLL_InitAlltoallV(MPI_Comm comm)
{
    //特别注意，comm通信子必须以节点数量为分配单位
    //目前拟支持每节点一个进程。进程内采用openmp。
    MPI_Comm_size(comm, &group_procn);
    MPI_Comm_rank(comm, &group_rank);
    Comm_group = comm;
    char *type = getenv("AlltoallvAlgorithmTYPE");
    if (strcmp(type, "LINEAR_SHIFT") == 0)
    {
        //基础节点间算法
        AlltoallvAlgorithm = LINEAR_SHIFT;
    }
    else if (strcmp(type, "XOR_EXCHANGE") == 0)
    {
        AlltoallvAlgorithm = XOR;
    }
    else if (strcmp(type, "OPT_BANDWIDTH") == 0)
    {
        AlltoallvAlgorithm = OPT_BANDWIDTH;
    }
    else if (strcmp(type, "LINEAR_SHIFT_REORDER") == 0)
    {
        AlltoallvAlgorithm = LINEAR_SHIFT_REORDER;
        char *Order_type = getenv("ORDER_TYPE");
        if (strcmp(Order_type, "reverse_order") == 0)
        {
            AlltoallvReorderType = reverse_order;
            //puts("check AlltoallvReorderType");
        }
        else if (strcmp(Order_type, "random_order") == 0)
        {
            //生成随机的步数顺序
            random_order_vec = (int *)malloc(sizeof(int) * (group_procn - 1)); //[1,group_procn-1]的重排数组
            void pjt_swap(int *a, int *b)
            {
                int c = *a;
                *a = *b;
                *b = c;
            }
            AlltoallvReorderType = random_order;
            random_order_vec = (int *)malloc(sizeof(int) * (group_procn - 1)); //[1,group_procn-1]的重排数组
            if (group_procn - 1 > 0)
            {
                for (int i = 0; i < group_procn - 1; i++)
                    random_order_vec[i] = i + 1;
                srand(time(NULL)); //让数变得随机
                if (group_rank == 0)
                {
                    for (int i = group_procn - 2; i >= 0; i--)
                    {
                        int b = rand() % (i + 1);
                        pjt_swap(&(random_order_vec[i]), &(random_order_vec[b]));
                    }
                }
                MPI_Bcast(random_order_vec, group_procn - 1, MPI_INT, 0, Comm_group);
            }
        }
        else if (strcmp(Order_type, "diavation_order") == 0)
        {
            //生成随机的步数顺序
            AlltoallvReorderType = diavation_order;
            int elem_size = (sizeof(float) + sizeof(int));
            diavation_order_vec = (char *)malloc(elem_size * (group_procn - 1));
            char swap_vec[elem_size];
            void pjt_swap(char *a, char *b)
            {
                memcpy(swap_vec, a, elem_size);
                memcpy(a, b, elem_size);
                memcpy(b, swap_vec, elem_size);
            }
            if (group_procn - 1 > 0)
            {
                for (int i = 0; i < group_procn - 1; i++)
                    (*(int *)(diavation_order_vec + i * elem_size + sizeof(float))) = i + 1;
                {
                    srand(time(NULL));
                    if (group_rank == 0)
                    {
                        for (int i = group_procn - 2; i >= 0; i--)
                        {
                            int b = rand() % (i + 1);
                            pjt_swap(diavation_order_vec + i * elem_size, diavation_order_vec + b * elem_size);
                        }
                    }
                    MPI_Bcast(diavation_order_vec, elem_size * (group_procn - 1), MPI_CHAR, 0, Comm_group);
                }
            }
            //初始顺序为随机序，随着alltoallv的多次调用，会在运行时改变顺序。
        }
        else if (strcmp(Order_type, "ring_distance_order") == 0)
        {
            AlltoallvReorderType = ring_distance_order;
        }
    }
    else if (strcmp(type, "RDMA_LINEAR_SHIFT") == 0)
    {
        //RDMA_LINEAR_SHIFT
        //从这里开始使用rdma
        AlltoallvAlgorithm = RDMA_LINEAR_SHIFT;
        rmt_disps = (uint64_t *)malloc(sizeof(uint64_t) * group_procn);
        {
            //生成固定的初始发送顺序
            AlltoallvReorderType = diavation_order;
            int elem_size = (sizeof(float) + sizeof(int));
            diavation_order_vec = (char *)malloc(elem_size * (group_procn - 1));
            char swap_vec[elem_size];
            void pjt_swap(char *a, char *b)
            {
                memcpy(swap_vec, a, elem_size);
                memcpy(a, b, elem_size);
                memcpy(b, swap_vec, elem_size);
            }
            if (group_procn - 1 > 0)
            {
                for (int i = 0; i < group_procn - 1; i++)
                    (*(int *)(diavation_order_vec + i * elem_size + sizeof(float))) = i + 1;
            }
        }
    }
    else if (strcmp(type, "RDMA_LINEAR_SHIFT_REORDER") == 0)
    {
        AlltoallvAlgorithm = RDMA_LINEAR_SHIFT_REORDER;
        rmt_disps = (uint64_t *)malloc(sizeof(uint64_t) * group_procn);
        char *Order_type = getenv("ORDER_TYPE");
        if (strcmp(Order_type, "random_order") == 0)
        {
            //生成随机的步数顺序
            AlltoallvReorderType = diavation_order;
            int elem_size = (sizeof(float) + sizeof(int));
            diavation_order_vec = (char *)malloc(elem_size * (group_procn - 1));
            char swap_vec[elem_size];
            void pjt_swap(char *a, char *b)
            {
                memcpy(swap_vec, a, elem_size);
                memcpy(a, b, elem_size);
                memcpy(b, swap_vec, elem_size);
            }
            if (group_procn - 1 > 0)
            {
                for (int i = 0; i < group_procn - 1; i++)
                    (*(int *)(diavation_order_vec + i * elem_size + sizeof(float))) = i + 1;
                // {
                //     srand(time(NULL));
                //     if (group_rank == 0)
                //     {
                //         for (int i = group_procn - 2; i >= 0; i--)
                //         {
                //             int b = rand() % (i + 1);
                //             pjt_swap(diavation_order_vec + i * elem_size, diavation_order_vec + b * elem_size);
                //         }
                //     }
                //     MPI_Bcast(diavation_order_vec, elem_size * (group_procn - 1), MPI_CHAR, 0, Comm_group);
                // }
            }
        }
    }
    else if (strcmp(type, "LINEAR_SHIFT_COMPRESS") == 0)
    {
        AlltoallvAlgorithm = LINEAR_SHIFT_COMPRESS;
    }
#ifdef PERFORMANCE_DEBUG
    statistics_time = (double *)malloc(sizeof(double) * global_procn);
    for (int i = 0; i < global_procn; i++)
        statistics_time[i] = 0.0;
#endif

    for (int i = 0; i < 2; i++)
    {
        AlltoallvBufSends[i] = (char *)(malloc(1 << 18));
        AlltoallvBufRecvs[i] = (char *)(malloc(1 << 18));
        ;
    }
    reqSV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    reqRV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    statusV = (MPI_Status *)malloc(sizeof(MPI_Status) * alltoallV_ongoingmsgN);
    total_time_MAX = (double *)malloc(sizeof(double) * group_procn);
    total_time_MIN = (double *)malloc(sizeof(double) * group_procn);
}

//分配一块通信内存
//目前只支持每节点放置一个进程
int GLECOLL_Register_SendRecv_Pair(void *sendbuf, uint64_t size1, void *recvbuf, uint64_t size2, MPI_Comm comm)
{
    glex_mem_handle_t local_mh;
    // puts("239");
    // printf("sendbuf = %p\n",sendbuf);
    // printf("ep = %p\n",_GLEXCOLL.ep);
    int ret = glex_register_mem(_GLEXCOLL.ep, sendbuf, size1,
                                GLEX_MEM_READ | GLEX_MEM_WRITE,
                                &local_mh);
                                // puts("242");

    MPI_Allgather((void *)&local_mh, sizeof(local_mh), MPI_CHAR,
                  send_mhs[iphmh_num], sizeof(local_mh), MPI_CHAR, comm);
    ret = glex_register_mem(_GLEXCOLL.ep, recvbuf, size2,
                            GLEX_MEM_READ | GLEX_MEM_WRITE,
                            &local_mh);
    MPI_Allgather((void *)&local_mh, sizeof(local_mh), MPI_CHAR,
                  recv_mhs[iphmh_num], sizeof(local_mh), MPI_CHAR, comm);
    iphmh_num++;
    MPI_Barrier(comm);
    return iphmh_num - 1;
}
void GLEXCOLL_Alltoallv_regist_recv_disps(uint64_t *rdispls)
{
    MPI_Alltoall(rdispls, 8, MPI_CHAR, rmt_disps, 8, MPI_CHAR, Comm_group);
    // printf("rank=%d rdispls[0]=%d rdispls[1]=%d\n",group_rank,rdispls[0],rdispls[1]);
    // printf("rank=%d rmt_disps[0]=%d rmt_disps[1]=%d\n",group_rank,rmt_disps[0],rmt_disps[1]);
}

void glexcoll_alltoallv_Linear_shift_REORDER_exchange_MPI(void *sendbuf,
                                                          uint64_t sendcounts[],
                                                          uint64_t sdispls[],
                                                          void *recvbuf,
                                                          uint64_t recvcounts[],
                                                          uint64_t rdispls[])
{
    static int call_count = 0;
    MPI_Request *reqSV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Request *reqRV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Status *statusV = (MPI_Status *)malloc(sizeof(MPI_Status) * alltoallV_ongoingmsgN);
    statistics_count = 0;
#ifdef PERFORMANCE_DEBUG
    double TIME_XE = MPI_Wtime();
#endif
    memcpy(recvbuf + rdispls[group_rank], sendbuf + sdispls[group_rank], sendcounts[group_rank]);

#ifdef PERFORMANCE_DEBUG
    statistics_time[statistics_count++] += (MPI_Wtime() - TIME_XE);
    //statistics_count++;
#endif
    //puts("start alltoallv");
    for (int i = 1; i < group_procn; i += alltoallV_ongoingmsgN)
    {
        int StartStep = 0, sign = 1;
        switch (AlltoallvReorderType)
        {
        case reverse_order:
            /* code */
            StartStep = group_procn - i;
            // if(group_rank ==0)
            //     printf("StartStep=%d\n",StartStep);
            sign = -1;
            break;
        case random_order:
            StartStep = random_order_vec[i - 1];
            sign = 1;
            break;
        case diavation_order:
            StartStep = *(int *)(diavation_order_vec + (i - 1) * 8 + 4);
            // if(group_rank == 0)
            //     printf("StartStep = %d\n",StartStep);
            sign = 1;
            break;
        case ring_distance_order:
            StartStep = group_procn - i;
            if ((i & 0x1) == 0)
            {
                //i为偶数
                StartStep = group_procn - i / 2;
            }
            else
            {
                //i为奇数
                StartStep = (i - 1) / 2 + 1;
            }
            // if(group_rank == 0)
            //     printf("StartStep = %d\n",StartStep);
            sign = 1;
            break;

        default:
            break;
        }
#ifdef PERFORMANCE_DEBUG
        double startT = MPI_Wtime();
#endif
        int tmp = 0;
        for (; tmp < alltoallV_ongoingmsgN; tmp++)
        {
            int Rs = StartStep + tmp * sign;
            if (Rs >= group_procn || Rs <= 0)
                break;
            int target = (group_rank + Rs) % group_procn;
            int source = (group_rank + group_procn - Rs) % group_procn;
            MPI_Isend(sendbuf + sdispls[target], sendcounts[target], MPI_CHAR, target, 0, Comm_group, reqSV + tmp);
            MPI_Irecv(recvbuf + rdispls[source], recvcounts[source], MPI_CHAR, source, 0, Comm_group, reqRV + tmp);
            //if(group_rank ==0)
        }
        if (tmp > 0)
        {
            MPI_Waitall(tmp, reqSV, statusV);
            MPI_Waitall(tmp, reqRV, statusV);
            //MPI_Barrier(MPI_COMM_WORLD);
        }
#ifdef PERFORMANCE_DEBUG
        double time_total = MPI_Wtime() - startT;
        statistics_time[statistics_count++] += time_total;
#endif
    }
    call_count++;
    if (call_count % 8 == 1 &&
        AlltoallvReorderType == diavation_order)
    {
        alltoallv_update_communication_order();
    }
    free(reqSV);
    free(reqRV);
    free(statusV);
}
//注意有一部分是char型，有一部分是整形
void glexcoll_alltoallv_Linear_shift_exchange_MPI(void *sendbuf,
                                                  uint64_t sendcounts[],
                                                  uint64_t sdispls[],
                                                  void *recvbuf,
                                                  uint64_t recvcounts[],
                                                  uint64_t rdispls[],
                                                  MPI_Datatype type)
{
    int elem_size = 0;
    MPI_Type_size(type, &elem_size);
    MPI_Request *reqSV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Request *reqRV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Status *statusV = (MPI_Status *)malloc(sizeof(MPI_Status) * alltoallV_ongoingmsgN);
    statistics_count = 0;
    // #ifdef PERFORMANCE_DEBUG
    //     double TIME_XE = MPI_Wtime();
    // #endif
    memcpy(recvbuf + elem_size * rdispls[group_rank],
           sendbuf + elem_size * sdispls[group_rank],
           elem_size * sendcounts[group_rank]);

    // #ifdef PERFORMANCE_DEBUG
    //     statistics_time[statistics_count++] += (MPI_Wtime() - TIME_XE);
    //     //statistics_count++;
    // #endif

    for (int shift = 1; shift < group_procn; shift += alltoallV_ongoingmsgN)
    {

        // #ifdef PERFORMANCE_DEBUG
        //         double startT = MPI_Wtime();
        // #endif
        int tmp = 0;
        for (; tmp < alltoallV_ongoingmsgN; tmp++)
        {
            int Rs = shift + tmp;
            if (Rs >= group_procn)
                break;
            int target = (group_rank + Rs) % group_procn;
            int source = (group_rank + group_procn - Rs) % group_procn;
            MPI_Isend(sendbuf + elem_size * sdispls[target],
                      sendcounts[target], type, target, 0, Comm_group, reqSV + tmp);
            MPI_Irecv(recvbuf + elem_size * rdispls[source],
                      recvcounts[source], type, source, 0, Comm_group, reqRV + tmp);
            //if(group_rank ==0)
        }
        if (tmp > 0)
        {
            MPI_Waitall(tmp, reqSV, statusV);
            MPI_Waitall(tmp, reqRV, statusV);
            //MPI_Barrier(MPI_COMM_WORLD);
        }
        // #ifdef PERFORMANCE_DEBUG
        //         double time_total = MPI_Wtime() - startT;
        //         // statistics_time[shift] += time_total;
        //         // statistics_count++;
        //         statistics_time[statistics_count++] += time_total;
        // #endif
    }
    free(reqSV);
    free(reqRV);
    free(statusV);
}
//一下方法适用于每节点多进程
//基本原理是通过整形数据压缩，将消息压缩成更小的消息
int alltoallv_single_data_size = 4;
void glexcoll_alltoallv_LINEAR_SHIFT_COMPRESS_MPI(void *sendbuf,
                                                  uint64_t sendcounts[],
                                                  uint64_t sdispls[],
                                                  void *recvbuf,
                                                  uint64_t recvcounts[],
                                                  uint64_t rdispls[])
{
//     MPI_Request *reqSV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
//     MPI_Request *reqRV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
//     MPI_Status *statusV = (MPI_Status *)malloc(sizeof(MPI_Status) * alltoallV_ongoingmsgN);
//     statistics_count = 0;
// #ifdef PERFORMANCE_DEBUG
//     double TIME_XE = MPI_Wtime();
// #endif
//     memcpy(recvbuf + alltoallv_single_data_size * rdispls[group_rank],
//            sendbuf + alltoallv_single_data_size * sdispls[group_rank],
//            alltoallv_single_data_size * sendcounts[group_rank]);

// #ifdef PERFORMANCE_DEBUG
//     statistics_time[statistics_count++] += (MPI_Wtime() - TIME_XE);
//     //statistics_count++;
// #endif
//     int source, target;
//     int prev_S, prev_T;
//     int shift = 1;
//     int bufid = shift & 1;
//     int prevbufid = bufid;
//     {
//         //对于shift=1时，把消息直接发送出去，不做任何压缩
//         {
//             target = (group_rank + shift) % group_procn;
//             source = (group_rank + group_procn - shift) % group_procn;
//             prev_S = source;
//             prev_T = target;
//             MPI_Isend(sendbuf + alltoallv_single_data_size * sdispls[target],
//                       sendcounts[target], MPI_INT, target, 0, Comm_group, reqSV);
//             MPI_Irecv(AlltoallvBufRecvs[bufid], recvcounts[source] + 10000, MPI_INT, source, 0, Comm_group, reqRV);
//         }
//     }
//     for (shift = 2; shift < group_procn; shift++)
//     {
//         int compressed_size = 0;
//         {
//             //压缩当前步消息
//             prevbufid = bufid;
//             bufid = shift & 1;
//             prev_S = source;
//             prev_T = target;
//             target = (group_rank + shift) % group_procn;
//             source = (group_rank + group_procn - shift) % group_procn;
//             compressed_size = p4nenc128v32((int *)(sendbuf + alltoallv_single_data_size * sdispls[target]),
//                                            sendcounts[target],
//                                            (int *)AlltoallvBufSends[bufid]);
//             compressed_size = sendcounts[target] * alltoallv_single_data_size;
//         }
//         {
//             MPI_Waitall(1, reqSV, statusV);
//             MPI_Waitall(1, reqRV, statusV);
//             // puts("check start");
//             MPI_Isend(AlltoallvBufSends[bufid], compressed_size, MPI_CHAR, target, 0, Comm_group, reqSV);
//             MPI_Irecv(AlltoallvBufRecvs[bufid], alltoallv_single_data_size * (recvcounts[source] + 10),
//                       MPI_CHAR, source, 0, Comm_group, reqRV);
//         }
//         // if (shift != 2)
//         // {
//         //     //将上一步消息解压缩到 recvbuf + alltoallv_single_data_size*rdispls[source]
//         //     p4ndec128v32((int *)AlltoallvBufRecvs[prevbufid],
//         //                  recvcounts[prev_S],
//         //                  recvbuf + alltoallv_single_data_size * rdispls[prev_S]);
//         // }
//         // else
//         // {
//         //     memcpy(recvbuf + alltoallv_single_data_size * rdispls[prev_S],
//         //            AlltoallvBufRecvs[prevbufid],
//         //            recvcounts[prev_S] * alltoallv_single_data_size);
//         // }
//     }
//     {
//         MPI_Waitall(1, reqSV, statusV);
//         MPI_Waitall(1, reqRV, statusV);
//         // puts("check");
//         // if (group_procn > 2)
//         //     p4ndec128v32((int *)AlltoallvBufRecvs[bufid],
//         //                  recvcounts[source],
//         //                  recvbuf + alltoallv_single_data_size * rdispls[source]);
//         // else
//         // {
//         //     memcpy(recvbuf + alltoallv_single_data_size * rdispls[source],
//         //            AlltoallvBufRecvs[bufid],
//         //            recvcounts[source] * alltoallv_single_data_size);
//         // }
//     }
//     free(reqSV);
//     free(reqRV);
//     free(statusV);
}

void glexcoll_alltoallv_XOR_exchange_MPI_PERFORMANCE_EXPLORE(void *sendbuf,
                                                             uint64_t sendcounts[],
                                                             uint64_t sdispls[],
                                                             void *recvbuf,
                                                             uint64_t recvcounts[],
                                                             uint64_t rdispls[])
{

    //for(int step_start = 1;step_start<group_procn;step_start+=alltoallV_ongoingmsgN)
    for (int step_start = 0; step_start < 26; step_start += alltoallV_ongoingmsgN)
    {

#ifdef PERFORMANCE_DEBUG
        double startT = MPI_Wtime();
#endif
        int tmp = 0;
        for (; tmp < alltoallV_ongoingmsgN; tmp++)
        {
            //if(step%2 ==0) step+=1;
            //if(step%2 ==1) step-=1;
            // if(step%2 ==1) step=31;
            // else step =32;
            // step = group_procn/2;
            int step = 1;
            int target = group_rank ^ step; //(group_rank + Rs)%group_procn;
            int source = target;            //(group_rank + group_procn - Rs)%group_procn;
            MPI_Isend(sendbuf + sdispls[target], sendcounts[target], MPI_CHAR, target, 0, Comm_group, reqSV + tmp);
            MPI_Irecv(recvbuf + rdispls[source], recvcounts[source], MPI_CHAR, source, 0, Comm_group, reqRV + tmp);
            //if(group_rank ==0)
        }

        if (tmp > 0)
        {
            MPI_Waitall(tmp, reqSV, statusV);
            MPI_Waitall(tmp, reqRV, statusV);
            //MPI_Barrier(MPI_COMM_WORLD);
        }

#ifdef PERFORMANCE_DEBUG
        double time_total = MPI_Wtime() - startT;
        statistics_time[step_start] += time_total;
        statistics_count++;
#endif
    }
}

void glexcoll_alltoallv_XOR_exchange_MPI(void *sendbuf,
                                         uint64_t sendcounts[],
                                         uint64_t sdispls[],
                                         void *recvbuf,
                                         uint64_t recvcounts[],
                                         uint64_t rdispls[])
{

    statistics_count = 0;

#ifdef PERFORMANCE_DEBUG
    double TIME_XE = MPI_Wtime();
#endif
    memcpy(recvbuf + rdispls[group_rank], sendbuf + sdispls[group_rank], sendcounts[group_rank]);

#ifdef PERFORMANCE_DEBUG
    statistics_time[0] += (MPI_Wtime() - TIME_XE);
    statistics_count++;
#endif

    for (int step_start = 1; step_start < group_procn; step_start += alltoallV_ongoingmsgN)
    {

#ifdef PERFORMANCE_DEBUG
        double startT = MPI_Wtime();
#endif
        int tmp = 0;
        for (; tmp < alltoallV_ongoingmsgN; tmp++)
        {
            int step = step_start + tmp;
            //if(step%2 ==0) step+=1;
            //if(step%2 ==1) step-=1;
            // if(step%2 ==1) step=31;
            // else step =32;
            // step = group_procn/2;
            if (step >= group_procn)
                break;
            int target = group_rank ^ step; //(group_rank + Rs)%group_procn;
            int source = target;            //(group_rank + group_procn - Rs)%group_procn;
            MPI_Isend(sendbuf + sdispls[target], sendcounts[target], MPI_CHAR, target, 0, Comm_group, reqSV + tmp);
            MPI_Irecv(recvbuf + rdispls[source], recvcounts[source], MPI_CHAR, source, 0, Comm_group, reqRV + tmp);
            //if(group_rank ==0)
        }

        if (tmp > 0)
        {
            MPI_Waitall(tmp, reqSV, statusV);
            MPI_Waitall(tmp, reqRV, statusV);
            //MPI_Barrier(MPI_COMM_WORLD);
        }

#ifdef PERFORMANCE_DEBUG
        double time_total = MPI_Wtime() - startT;
        statistics_time[step_start] += time_total;
        statistics_count++;
#endif
    }
}

void glexcoll_alltoallv_OPT_BINDWIDTH_exchange_MPI(void *sendbuf,
                                                   uint64_t sendcounts[],
                                                   uint64_t sdispls[],
                                                   void *recvbuf,
                                                   uint64_t recvcounts[],
                                                   uint64_t rdispls[])
{
    MPI_Request *reqSV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Request *reqRV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Status *statusV = (MPI_Status *)malloc(sizeof(MPI_Status) * alltoallV_ongoingmsgN);
    //printf("group_procn=%d\n",group_procn & 0x1);
    if ((group_procn & 0x1) == 0)
    {
        //偶数进程情况
        //puts("check position");
        for (int step_start = 0; step_start < group_procn; step_start += alltoallV_ongoingmsgN)
        {
#ifdef PERFORMANCE_DEBUG
            double startT = MPI_Wtime();
#endif
            int tmp = 0;
            for (; tmp < alltoallV_ongoingmsgN; tmp++)
            {
                int step = step_start + tmp;
                if (step >= group_procn)
                    break;
                int target;
                if (group_rank % 2 == 0)
                {
                    target = ((group_rank >> 1) + step) % group_procn;
                }
                else
                {
                    target = ((group_procn >> 1) + ((group_rank - 1) >> 1) + step) % (group_procn);
                }
                int source;
                int i1 = (group_rank - step + group_procn) % group_procn;
                if (i1 < (group_procn >> 1))
                {
                    source = (2 * i1) % group_procn;
                }
                else
                {
                    source = ((i1 - group_procn / 2) * 2 + 1) % group_procn;
                }
                //printf("step=%d rank=%d source=%d target=%d\n",step,group_rank,source,target);
                MPI_Isend(sendbuf + sdispls[target], sendcounts[target], MPI_CHAR, target, flag_control, Comm_group, reqSV + tmp);
                MPI_Irecv(recvbuf + rdispls[source], recvcounts[source], MPI_CHAR, source, flag_control++, Comm_group, reqRV + tmp);
                //if(group_rank ==0)
            }
            if (tmp > 0)
            {
                MPI_Waitall(tmp, reqSV, statusV);
                MPI_Waitall(tmp, reqRV, statusV);
                //MPI_Barrier(MPI_COMM_WORLD);
                // if(group_rank == 0)
                //     printf("check%d\n",step_start);
            }
#ifdef PERFORMANCE_DEBUG
            double time_total = MPI_Wtime() - startT;
            statistics_time[step_start] += time_total;
            statistics_count++;
#endif
        }
    }
    else
    {
        //奇数进程情况
    }

    free(reqSV);
    free(reqRV);
    free(statusV);
}

//一次带有性能统计的alltoallv传输。通过该传输可以生成Linear Shift的防拥塞调度
void GLEXCOLL_Alltoallv_Performance_detection(void *sendbuf,
                                              uint64_t sendcounts[],
                                              uint64_t sdispls[],
                                              void *recvbuf,
                                              uint64_t recvcounts[],
                                              uint64_t rdispls[])
{
    static int call_count = 0;
    MPI_Request *reqSV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Request *reqRV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoallV_ongoingmsgN);
    MPI_Status *statusV = (MPI_Status *)malloc(sizeof(MPI_Status) * alltoallV_ongoingmsgN);
    statistics_count = 0;
#ifdef PERFORMANCE_DEBUG
    double TIME_XE = MPI_Wtime();
#endif
    memcpy(recvbuf + rdispls[group_rank], sendbuf + sdispls[group_rank], sendcounts[group_rank]);

#ifdef PERFORMANCE_DEBUG
    statistics_time[statistics_count++] += (MPI_Wtime() - TIME_XE);
    //statistics_count++;
#endif
    //puts("start alltoallv");
    for (int i = 1; i < group_procn; i += alltoallV_ongoingmsgN)
    {
        int StartStep = 0, sign = 1;
        StartStep = *(int *)(diavation_order_vec + (i - 1) * 8 + 4);
        // if(group_rank == 0)
        //     printf("StartStep = %d\n",StartStep);
        sign = 1;
#ifdef PERFORMANCE_DEBUG
        double startT = MPI_Wtime();
#endif
        int tmp = 0;
        for (; tmp < alltoallV_ongoingmsgN; tmp++)
        {
            int Rs = StartStep + tmp * sign;
            if (Rs >= group_procn || Rs <= 0)
                break;
            int target = (group_rank + Rs) % group_procn;
            int source = (group_rank + group_procn - Rs) % group_procn;
            MPI_Isend(sendbuf + sdispls[target], sendcounts[target], MPI_CHAR, target, 0, Comm_group, reqSV + tmp);
            MPI_Irecv(recvbuf + rdispls[source], recvcounts[source], MPI_CHAR, source, 0, Comm_group, reqRV + tmp);
            //if(group_rank ==0)
        }
        if (tmp > 0)
        {
            MPI_Waitall(tmp, reqSV, statusV);
            MPI_Waitall(tmp, reqRV, statusV);
            //MPI_Barrier(MPI_COMM_WORLD);
        }
#ifdef PERFORMANCE_DEBUG
        double time_total = MPI_Wtime() - startT;
        statistics_time[statistics_count++] += time_total;
#endif
    }
    {
        alltoallv_update_communication_order();
    }
    free(reqSV);
    free(reqRV);
    free(statusV);
}
void GLEXCOLL_Alltoallv_RDMA_REORDER(int mh_pair, void *sendbuf, void *recvbuf, uint64_t sendcounts[], uint64_t sdispls[])
{
    //puts("check start");
    static struct glex_rdma_req rdma_req;
    static struct glex_rdma_req *bad_rdma_req;

    int on_going_count = 0;

    for (int i = 1; i < group_procn; ++i)
    {
        int StartStep = *(int *)(diavation_order_vec + (i - 1) * 8 + 4);
        int target = (group_rank + StartStep) % group_procn;
        // if(group_rank == 1)
        //     printf("GLEXCOLL_Alltoallv_RDMA_REORDER: target=%d\n",target);
        {
            rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
            // //puts("check 225");
            rdma_req.local_mh.v = send_mhs[mh_pair][group_rank].v;
            rdma_req.local_offset = sdispls[target];
            rdma_req.len = sendcounts[target];
            rdma_req.rmt_mh.v = recv_mhs[mh_pair][target].v;
            rdma_req.rmt_offset = rmt_disps[target];
            // //puts("check 231");
            rdma_req.type = GLEX_RDMA_TYPE_PUT;
            rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
            rdma_req.rmt_evt.cookie[1] = 999;
            // rdma_req.local_evt.cookie[0] = 998;
            // rdma_req.local_evt.cookie[1] = 997;
            rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
            rdma_req.flag = GLEX_FLAG_REMOTE_EVT; //| GLEX_FLAG_LOCAL_EVT;
            if (i % alltoallV_ongoingmsgN == 0)
                rdma_req.flag |= GLEX_FLAG_FENCE;
            rdma_req.next = NULL;
            int ret;
            while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
            {
            }
            if (ret != GLEX_SUCCESS)
            {
                if (ret == GLEX_INVALID_PARAM)
                    printf("%d, _rdma() 非法参数", global_rank);
                printf("_rdma(), return: %d\n", ret);
                exit(1);
            }
        }
    }
    memcpy(recvbuf + rmt_disps[group_rank], sendbuf + sdispls[group_rank], sendcounts[group_rank]);
    for (int i = 1; i < group_procn; ++i)
    {
        static glex_event_t *event;
        while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            ;
        if (event->cookie[1] != 999)
        {
            puts("check cookie error");
            exit(0);
        }
    }
    glex_discard_probed_event(_GLEXCOLL.ep);
    //printf("rank=%d rmt_disps[group_rank]=%d sdispls[group_rank]=%d sendcounts[group_rank]=%d\n",
    // group_rank,
    // rmt_disps[group_rank],
    // sdispls[group_rank],
    // sendcounts[group_rank]);
}
void GLEXCOLL_Alltollv_RDMA(int mh_pair, void *sendbuf, void *recvbuf, uint64_t sendcounts[], uint64_t sdispls[])
{
    //puts("check start");
    static struct glex_rdma_req rdma_req;
    static struct glex_rdma_req *bad_rdma_req;

    int on_going_count = 0;
    for (int i = 1; i < group_procn; ++i)
    {
        int target = (group_rank + i) % group_procn;
        {
            rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
            // //puts("check 225");
            rdma_req.local_mh.v = send_mhs[mh_pair][group_rank].v;
            rdma_req.local_offset = sdispls[target];
            rdma_req.len = sendcounts[target];
            rdma_req.rmt_mh.v = recv_mhs[mh_pair][target].v;
            rdma_req.rmt_offset = rmt_disps[target];
            // //puts("check 231");
            rdma_req.type = GLEX_RDMA_TYPE_PUT;
            rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
            rdma_req.rmt_evt.cookie[1] = 999;
            // rdma_req.local_evt.cookie[0] = 998;
            // rdma_req.local_evt.cookie[1] = 997;
            rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
            rdma_req.flag = GLEX_FLAG_REMOTE_EVT; //| GLEX_FLAG_LOCAL_EVT;
            if (i % 2 == 0)
                rdma_req.flag |= GLEX_FLAG_FENCE;
            rdma_req.next = NULL;
            int ret;
            while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
            {
            }
            if (ret != GLEX_SUCCESS)
            {
                if (ret == GLEX_INVALID_PARAM)
                    printf("%d, _rdma() 非法参数", global_rank);
                printf("_rdma(), return: %d\n", ret);
                exit(1);
            }
        }
    }
    memcpy(recvbuf + rmt_disps[group_rank], sendbuf + sdispls[group_rank], sendcounts[group_rank]);
    for (int i = 1; i < group_procn; ++i)
    {
        static glex_event_t *event;
        while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            ;
        if (event->cookie[1] != 999)
        {
            puts("check cookie error");
            exit(0);
        }
    }
    glex_discard_probed_event(_GLEXCOLL.ep);
    //printf("rank=%d rmt_disps[group_rank]=%d sdispls[group_rank]=%d sendcounts[group_rank]=%d\n",
    // group_rank,
    // rmt_disps[group_rank],
    // sdispls[group_rank],
    // sendcounts[group_rank]);
}

void alltoallv_update_communication_order()
{
    int pjt_cmp(const void *a, const void *b)
    {
        // if((*(float *)a) > (*(float *)b);)//从大到小
        //     return 1;
        if ((*(float *)b) > (*(float *)a)) //从小到大
            return 1;
        return 0;
    }
    MPI_Barrier(Comm_group);
    MPI_Request reqa, reqb;
    MPI_Status stat;
    // {
    //     MPI_Iallreduce(statistics_time, total_time_MAX, statistics_count, MPI_DOUBLE, MPI_MAX, Comm_group, &reqa);
    //     MPI_Iallreduce(statistics_time, total_time_MIN, statistics_count, MPI_DOUBLE, MPI_MIN, Comm_group, &reqb);
    //     MPI_Wait(&reqa, &stat);
    //     MPI_Wait(&reqb, &stat);
    //     for (int i = 1; i < statistics_count; i++)
    //     {
    //         float diavation = total_time_MAX[i] - total_time_MIN[i];
    //         *(float *)(diavation_order_vec + (i - 1) * DIAVATION_VEC_ELEM_SIZE) = diavation;
    //         statistics_time[i] = 0.0;
    //     }
    // }
    {
        MPI_Allreduce(statistics_time, total_time_MAX, statistics_count, MPI_DOUBLE, MPI_SUM, Comm_group);
        for (int i = 1; i < statistics_count; i++)
        {
            float diavation = total_time_MAX[i] / statistics_count;
            *(float *)(diavation_order_vec + (i - 1) * DIAVATION_VEC_ELEM_SIZE) = diavation;
            statistics_time[i] = 0.0;
        }
    }

    statistics_time[0] = 0.0;
    qsort(diavation_order_vec, statistics_count - 1, DIAVATION_VEC_ELEM_SIZE, &pjt_cmp);

    // if(group_rank == 0)
    //     puts("sort check");
}
void GLEXCOLL_Alltoallv(void *sendbuf,
                        uint64_t sendcounts[],
                        uint64_t sdispls[],
                        void *recvbuf,
                        uint64_t recvcounts[],
                        uint64_t rdispls[],
                        MPI_Datatype type)
{
    switch (AlltoallvAlgorithm)
    {
    case LINEAR_SHIFT:
        //puts("LINEAR_SHIFT");
        glexcoll_alltoallv_Linear_shift_exchange_MPI(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, type);
        break;
    case XOR:
        //puts("XOR");
        glexcoll_alltoallv_XOR_exchange_MPI(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls);
        break;
    case OPT_BANDWIDTH:
        //puts("OPT_BANDWIDTH");
        glexcoll_alltoallv_OPT_BINDWIDTH_exchange_MPI(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls);
        break;
    case LINEAR_SHIFT_REORDER:
        glexcoll_alltoallv_Linear_shift_REORDER_exchange_MPI(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls);
        break;
    case LINEAR_SHIFT_COMPRESS:
        glexcoll_alltoallv_LINEAR_SHIFT_COMPRESS_MPI(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls);
    default:
        break;
    }
}

void write_out_print_alltoallv_performance_data(int loopn, double *Step_max_sum)
{
    //loopn=1;
    //statistics_count/=loopn;
    double *total_time = (double *)malloc(sizeof(double) * statistics_count);
    double *all_statistics_data = 0;

    if (group_rank == 0)
        all_statistics_data = (double *)malloc(sizeof(double) * group_procn * statistics_count);
    MPI_Reduce(statistics_time, total_time, statistics_count, MPI_DOUBLE, MPI_SUM, 0, Comm_group);
    MPI_Reduce(statistics_time, total_time_MAX, statistics_count, MPI_DOUBLE, MPI_MAX, 0, Comm_group);
    MPI_Reduce(statistics_time, total_time_MIN, statistics_count, MPI_DOUBLE, MPI_MIN, 0, Comm_group);
    MPI_Gather(statistics_time, statistics_count, MPI_DOUBLE, all_statistics_data, statistics_count, MPI_DOUBLE, 0, Comm_group);
    MPI_Barrier(Comm_group);
    (*Step_max_sum) = 0.0;
    for (int i = 0; i < statistics_count; i++)
    {
        total_time[i] /= (loopn * group_procn);
        total_time_MAX[i] /= loopn;
        total_time_MIN[i] /= loopn;
        (*Step_max_sum) += total_time_MAX[i];
    }
    if (group_rank == 0)
        for (int i = 0; i < group_procn; i++)
        {
            for (int j = 0; j < statistics_count; j++)
                all_statistics_data[i * statistics_count + j] /= loopn;
        }
    //puts("XXXXXXXXXXX");
    {
        char outputF[60] = "alltoallv_performance.txt";
        if (group_rank == 0)
        {
            FILE *fp = fopen(outputF, "ab+");
            fprintf(fp, "%d\n", statistics_count);
            for (int i = 0; i < statistics_count; i++)
                fprintf(fp, "%lf\t%lf\t%lf \n", 1e6 * total_time_MIN[i], 1e6 * total_time[i], 1e6 * total_time_MAX[i]);
            fputs("\n", fp);
            fclose(fp);
        }
    }
    {
        char outputF[60] = "alltoallv_all_statistics.txt";
        if (group_rank == 0)
        {
            FILE *fp = fopen(outputF, "ab+");
            fprintf(fp, "%d %d\n", group_procn, statistics_count);
            for (int i = 0; i < group_procn; i++)
            {
                for (int j = 0; j < statistics_count; j++)
                {
                    fprintf(fp, "%lf ", 1e6 * all_statistics_data[i * statistics_count + j]);
                }
                fputs("\n", fp);
            }
            fputs("\n", fp);
            fclose(fp);
        }
    }

    free(total_time);
    MPI_Barrier(Comm_group);
}
void GLEXCOLL_AlltoallVFinalize()
{

#ifdef PERFORMANCE_DEBUG
    free(statistics_time);
#endif
    free(reqSV);
    free(reqRV);
    free(statusV);
    free(total_time_MAX);
    free(total_time_MIN);
}