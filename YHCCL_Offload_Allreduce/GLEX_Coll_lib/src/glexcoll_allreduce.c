#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <endian.h>
#include <sched.h>
#include <errno.h>
#include <mpi.h>
#include <stdint.h>
#include <omp.h>
#include <sys/mman.h>
#include <pthread.h>
#include "glex.h"
#include "glexcoll.h"
#include "glexallreduce.h"

static int allreduce_rank;
static int allreduce_procn;
static int allreduce_intra_rank;
static int allreduce_intra_procn;
static int allreduce_inter_rank;
static int allreduce_inter_procn;

static MPI_Comm allreduce_comm;
int allreduce_ongoin_msgN = 1;
int allreduce_slice_num = (1 << 10);

int SoftWare_Allreduce_Algorithm = Self_tuned_two_dimentional;

#define ALLREDUCE_TAG 996
void recursive_doubling_double_sum(double *sendbuf, double *recvbuf, int num)
{
    // puts("check pjt");
    int size = num * sizeof(double);
    static MPI_Request req_send;
    static MPI_Request req_recv;
    static MPI_Status status;
    double *tmpbuf = (double *)allreduce_sendbuf;
    memcpy(recvbuf, sendbuf, size);

    int stepN = allreduce_recursive_doubling_stepn;
    int step = 1;
    for (int i = 0; i < stepN; i++)
    {
        int target = (allreduce_rank ^ step);
        if (target < allreduce_procn)
        {
            MPI_Isend(recvbuf, size, MPI_CHAR, target, ALLREDUCE_TAG, allreduce_comm, &req_send);
            MPI_Irecv(tmpbuf, size, MPI_CHAR, target, ALLREDUCE_TAG, allreduce_comm, &req_recv);
            MPI_Wait(&req_recv, &status);
            MPI_Wait(&req_send, &status);
            // {
            //     //做一次RDMA发送和接收
            //     memccpy(allreduce_sendbuf,recv)
            // }
            // #pragma omp parallel for
            for (int j = 0; j < num; j++)
            {
                recvbuf[j] += tmpbuf[j];
            }
        }
        step *= 2;
    }
}

void PJT_discard_event(glex_ep_handle_t ep)
{
    _GLEXCOLL.event_credit -= 1; //= Event_CREDIT_MAX;
    if (_GLEXCOLL.event_credit <= 0)
    {
        glex_discard_probed_event(ep);
        _GLEXCOLL.event_credit = Event_CREDIT_MAX;
    }
}

static int allreduce_event_vec[100];
static int allreduce_event_count;
extern int allreduce_buf_length;
void recursive_doubling_double_sum_rdma(double *sendbuf, int startShift, double *recvbuf, int num)
{
    // static int round = 0;
    int start_offset = startShift;
    // round *allreduce_recursive_doubling_stepn *num * sizeof(double) % (1 << 11);
    //  start_offset=0;
    allreduce_event_count++;
    // puts("check");
    static struct glex_rdma_req rdma_req;
    static struct glex_rdma_req *bad_rdma_req;
    static glex_event_t *event;
    glex_ret_t ret;

    // puts("check pjt");
    int size = num * sizeof(double);
    static MPI_Request req_send;
    static MPI_Request req_recv;
    static MPI_Status status;
    double *tmpbuf = (double *)allreduce_sendbuf;
    if (intra_procn == 1)
    {
        memcpy(((void *)allreduce_sendbuf) + startShift, sendbuf, size);
    }

    int stepN = allreduce_recursive_doubling_stepn;
    int step = 1;

    double *sendb;

    for (int i = 0; i < stepN; i++)
    {
        int target = (allreduce_inter_rank ^ step);
        // printf("%d's target = %d\n",allreduce_inter_rank,target);
        if (target < allreduce_inter_procn)
        {
            // MPI_Isend(recvbuf, size, MPI_CHAR, target, ALLREDUCE_TAG, allreduce_comm, &req_send);
            // MPI_Irecv(tmpbuf, size, MPI_CHAR, target, ALLREDUCE_TAG, allreduce_comm, &req_recv);
            // MPI_Wait(&req_recv, &status);
            // MPI_Wait(&req_send, &status);
            {
                //做一次RDMA发送和接收
                // memccpy(allreduce_sendbuf,recv)
                {
                    //第一步把消息发给对方
                    {
                        rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                        rdma_req.local_mh.v = send_mhs[allreduce_send_recv_pair][allreduce_inter_rank].v;
                        rdma_req.local_offset = start_offset + i * size;
                        rdma_req.len = num * sizeof(double);
                        rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
                        rdma_req.rmt_offset = start_offset + i * size;
                        // printf("_ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n",_ReduceTreeS[_TreeID].DownReduceRank - 1);
                        rdma_req.type = GLEX_RDMA_TYPE_PUT;
                        rdma_req.rmt_evt.cookie[0] = 997;
                        rdma_req.rmt_evt.cookie[1] = i + 1;
                        // rdma_req.local_evt.cookie[0] = 996;
                        // rdma_req.local_evt.cookie[1] = i + 1;
                        rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                        rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
                        if (i == 0)
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
                    //第二步等待发送方的事件，注意这里有可能实际等到的是发送方下一次allreduce传输过来的消息，此时通过start_offset来区分。
                    while (allreduce_event_vec[i] < allreduce_event_count)
                    {
                        while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                            ;
                        allreduce_event_vec[event->cookie[1] - 1] += 1;
                        PJT_discard_event(_GLEXCOLL.ep);
                    }
                    //第三步规约消息
                    double *tmpbuf = allreduce_recvbuf + start_offset + i * (num * sizeof(double));
                    // printf("rank = %d *tmpbuf = %f asend_buf=%f\n",allreduce_rank,*tmpbuf,*(double *)allreduce_sendbuf);
                    sendb = allreduce_sendbuf + start_offset + (i + 1) * (num * sizeof(double));
                    double *A = allreduce_sendbuf + start_offset + (i) * (num * sizeof(double));
// #pragma omp parallel for num_threads(4)
#pragma omp simd
                    for (int j = 0; j < num; j++)
                    {
                        sendb[j] = A[j] + tmpbuf[j];
                    }
                }
            }
        }
        step *= 2;
    }
    memcpy(recvbuf, sendb, size);
}
#define compiler_barrier() __asm__ __volatile__("" \
                                                :  \
                                                :  \
                                                : "memory")
#define lfence() __asm__ __volatile__("lfence" \
                                      :        \
                                      :        \
                                      : "memory")

#define mfence() __asm__ __volatile__("mfence" \
                                      :        \
                                      :        \
                                      : "memory")
double *recursive_doubling_double_sum_rdma_avoidEvent(int startshift, int num)
{
    static int round = 996;
    //消息直接从alltoallsendbuf到recvbug
    // start_offset=0;
    round++;
    // puts("check");
    static struct glex_rdma_req rdma_req;
    static struct glex_rdma_req *bad_rdma_req;
    static glex_event_t *event;
    glex_ret_t ret;

    // puts("check pjt");
    int size0 = num * sizeof(double);
    static MPI_Request req_send;
    static MPI_Request req_recv;
    static MPI_Status status;

    int stepN = allreduce_recursive_doubling_stepn;
    int size = size0 + 4;
    for (int i = 0; i < stepN; i++)
    {
        *(volatile int *)(allreduce_sendbuf + startshift + size0 + i * size) = round;
    }
    int step = 1;

    // double *sendb;
    double *C;
    for (int i = 0; i < stepN; i++)
    {

        int target = (allreduce_inter_rank ^ step);
        //     // printf("%d's target = %d\n",allreduce_inter_rank,target);
        if (target < allreduce_inter_procn)
        {
            // MPI_Isend(recvbuf, size, MPI_CHAR, target, ALLREDUCE_TAG, allreduce_comm, &req_send);
            // MPI_Irecv(tmpbuf, size, MPI_CHAR, target, ALLREDUCE_TAG, allreduce_comm, &req_recv);
            // MPI_Wait(&req_recv, &status);
            // MPI_Wait(&req_send, &status);
            //         {
            //             //做一次RDMA发送和接收
            //             // memccpy(allreduce_sendbuf,recv)
            //             {
            //                 //第一步把消息发给对方
            {
                rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                rdma_req.local_mh.v = send_mhs[allreduce_send_recv_pair][allreduce_inter_rank].v;
                rdma_req.local_offset = startshift + i * size;
                rdma_req.len = size;
                rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
                rdma_req.rmt_evt.cookie[1] = i + 1;
                rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
                rdma_req.rmt_offset = startshift + i * size;
                // printf("_ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n",_ReduceTreeS[_TreeID].DownReduceRank - 1);
                rdma_req.type = GLEX_RDMA_TYPE_PUT;
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                rdma_req.flag = GLEX_FLAG_REMOTE_EVT; // |GLEX_FLAG_LOCAL_EVT;
                if (i == 0)
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
            //第二步等待target发送过来的消息
            volatile int *waitp = allreduce_recvbuf + startshift + (i + 1) * size - 4;

            // while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            //     ;
            while (*waitp != round)
                ;
            __sync_synchronize();
            // printf("%d check recv %d round:%d\n", inter_rank, *waitp, round);
            //                 while (allreduce_event_vec[i] < round)
            //                 {
            //                     while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            //                         ;
            //                     allreduce_event_vec[event->cookie[1] - 1] += 1;
            //                     PJT_discard_event(_GLEXCOLL.ep);
            //                 }
            //第三步规约消息
            // A为alltoallrecvbuf的0号部分，B为为alltoallrecvbuf的i号部分
            double *A = allreduce_sendbuf + startshift + i * size;
            double *B = allreduce_recvbuf + startshift + i * size;
            C = allreduce_sendbuf + startshift + (i + 1) * size;
            {
#pragma omp simd
                for (int j = 0; j < num; j++)
                {
                    C[j] = A[j] + B[j];
                    // if (C[j] < 63.999 || C[j] > 64.0001)
                    // {
                    //     printf("rank %d val:C[j] %f A[j] %f B[j] %f\n", inter_rank, C[j], A[j], B[j]);
                    //     printf("step %d error ", round);
                    // }
                }
            }
            // *(int *)(C + size0) = round;
            //                 double *tmpbuf = allreduce_recvbuf + start_offset + i * (num * sizeof(double));
            //                 // printf("rank = %d *tmpbuf = %f asend_buf=%f\n",allreduce_rank,*tmpbuf,*(double *)allreduce_sendbuf);
            //                 sendb = allreduce_sendbuf + (i + 1) * (num * sizeof(double));
            //                 double *A = allreduce_sendbuf + (i) * (num * sizeof(double));
            //                 for (int j = 0; j < num; j++)
            //                 {
            //                     sendb[j] = A[j] + tmpbuf[j];
            //                 }
            //             }
            //         }
            //         // #pragma omp parallel for
        }
        step *= 2;
    }
    return C;
    // memcpy(recvbuf, allreduce_sendbuf, num * sizeof(double));
}
double *recursive_doubling_double_sum_imm_rdma(int startshift, int num)
{
    static int round = 1;
    //消息直接从alltoallsendbuf到recvbug
    // start_offset=0;
    round++;
    // puts("check");
    static struct glex_imm_rdma_req rdma_req;
    static struct glex_imm_rdma_req *bad_rdma_req;
    glex_ret_t ret;

    // puts("check pjt");
    int size0 = num * sizeof(double);

    int stepN = allreduce_recursive_doubling_stepn;
    int size = size0 + 4;
    for (int i = 0; i < stepN; i++)
        *(volatile int *)(allreduce_sendbuf + startshift + size0 + i * size) = round;
    int step = 1;

    // double *sendb;
    volatile double *C;
    for (int i = 0; i < stepN; i++)
    {

        int target = (allreduce_inter_rank ^ step);
        //     // printf("%d's target = %d\n",allreduce_inter_rank,target);
        if (target < allreduce_inter_procn)
        {
            //                 //第一步把消息发给对方
            {
                rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                rdma_req.data = allreduce_sendbuf + startshift + i * size;
                rdma_req.len = size;
                rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
                rdma_req.rmt_offset = startshift + i * size;
                // printf("_ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n",_ReduceTreeS[_TreeID].DownReduceRank - 1);
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                rdma_req.flag = NULL;
                if (i == 0)
                    rdma_req.flag |= GLEX_FLAG_FENCE;
                rdma_req.next = NULL;
                int ret;
                while ((ret = glex_imm_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
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
            //第二步等待target发送过来的消息
            volatile int *waitp = allreduce_recvbuf + startshift + (i + 1) * size - 4;

            // while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            //     ;
            while (*waitp != round)
                ;
            __sync_synchronize();

            // printf("%d check recv %d round:%d\n", inter_rank, *waitp, round);
            //                 while (allreduce_event_vec[i] < round)
            //                 {
            //                     while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            //                         ;
            //                     allreduce_event_vec[event->cookie[1] - 1] += 1;
            //                     PJT_discard_event(_GLEXCOLL.ep);
            //                 }
            //第三步规约消息
            // A为alltoallrecvbuf的0号部分，B为为alltoallrecvbuf的i号部分
            volatile double *A = allreduce_sendbuf + startshift + i * size;
            volatile double *B = allreduce_recvbuf + startshift + i * size;
            C = allreduce_sendbuf + startshift + (i + 1) * size;
            {
                for (int j = 0; j < num; j++)
                {
                    C[j] = A[j] + B[j];
                    // if (C[j] < 63.999 || C[j] > 64.0001)
                    // {
                    //     printf("rank %d val:C[j] %f A[j] %f B[j] %f\n", inter_rank, C[j], A[j], B[j]);
                    //     printf("step %d error ", round);
                    // }
                }
            }
            // *(int *)(C + size0) = round;
            //                 double *tmpbuf = allreduce_recvbuf + start_offset + i * (num * sizeof(double));
            //                 // printf("rank = %d *tmpbuf = %f asend_buf=%f\n",allreduce_rank,*tmpbuf,*(double *)allreduce_sendbuf);
            //                 sendb = allreduce_sendbuf + (i + 1) * (num * sizeof(double));
            //                 double *A = allreduce_sendbuf + (i) * (num * sizeof(double));
            //                 for (int j = 0; j < num; j++)
            //                 {
            //                     sendb[j] = A[j] + tmpbuf[j];
            //                 }
            //             }
            //         }
            //         // #pragma omp parallel for
        }
        step *= 2;
    }
    return C;
    // memcpy(recvbuf, allreduce_sendbuf, num * sizeof(double));
}

extern int mmin(int a, int b);
void recursive_doubling_double_slicing_sum(double *sendbuf, double *recvbuf, int num)
{
    //     //puts("check pjt");
    //     int size = num*sizeof(double);
    //     static MPI_Request req_send[32];
    //     static MPI_Request req_recv[32];
    //     static MPI_Status  status[32];
    //     double * tmpbuf = (double *)malloc(size);
    //     memcpy(recvbuf,sendbuf,size);

    //     int stepN=0;
    //     int step=allreduce_procn;
    //     while (step > 1){
    //         stepN++;
    //         step/=2;
    //     }
    //     for(int i = 0;i<stepN;i++)
    //     {
    //         int BigBlockNum = allreduce_slice_num*allreduce_ongoin_msgN;
    //         int target = (allreduce_rank ^ step);
    //         for(int startshift = 0;startshift < num;startshift+=BigBlockNum)
    //         {
    //             extern int mmin(int a,int b);
    //             int shift_end=mmin(num,startshift+BigBlockNum);
    //             int s=0,shift=startshift;
    //             while (shift < shift_end)
    //             {
    //                 int local_num = mmin(allreduce_slice_num,shift_end - shift);
    //                 MPI_Isend(recvbuf+shift,local_num,MPI_DOUBLE,target,s,allreduce_comm,&(req_send[s]));
    //                 MPI_Irecv(tmpbuf+shift,local_num,MPI_DOUBLE,target,s,allreduce_comm,&(req_recv[s]));
    //                 /* code */
    //                 shift+=allreduce_slice_num;
    //                 s++;
    //             }
    //             MPI_Waitall(s,req_send,status);
    //             MPI_Waitall(s,req_recv,status);
    //         }
    //         //int count_tmp = mmin(num,BigBlockNum);
    //         //接下来进行第一波消息发�
    // #pragma omp parallel for
    //         for(int j = 0;j<num;j++)
    //         {
    //             recvbuf[j]+=tmpbuf[j];
    //         }
    //         step*=2;
    //     }
    //     free(tmpbuf);

    // puts("check pjt");
    int size = num * sizeof(double);
    static MPI_Request req_send;
    static MPI_Request req_recv;
    static MPI_Status status;
    double *tmpbuf = (double *)malloc(size);
    memcpy(recvbuf, sendbuf, size);

    int stepN = 0;
    int step = allreduce_procn;
    while (step > 1)
    {
        stepN++;
        step /= 2;
    }
    // printf("%d\n",stepN);
    for (int i = 0; i < stepN; i++)
    {
        int target = (allreduce_rank ^ step);
        // printf("%d %d\n",allreduce_rank,target);
        int last_start = 0, start = 0;
        int local_size = mmin(num, allreduce_slice_num);
        int last_local_size = local_size;
        MPI_Isend(recvbuf, local_size, MPI_DOUBLE, target, start, allreduce_comm, &req_send);
        MPI_Irecv(tmpbuf, local_size, MPI_DOUBLE, target, start, allreduce_comm, &req_recv);
        // for(int start=0;start<num;start+=allreduce_slice_num)
        start += local_size;
        while (start < num)
        {
            // puts("x");
            MPI_Wait(&req_recv, &status);
            MPI_Wait(&req_send, &status);

            last_local_size = local_size;
            local_size = mmin(num - start, allreduce_slice_num);
            MPI_Isend(&(recvbuf[start]), local_size, MPI_DOUBLE, target, start, allreduce_comm, &req_send);
            MPI_Irecv(&(tmpbuf[start]), local_size, MPI_DOUBLE, target, start, allreduce_comm, &req_recv);
            //开始上一轮完成传输的数据的规�
            {
#pragma omp parallel for
                for (int j = 0; j < last_local_size; j++)
                {
                    recvbuf[last_start + j] += tmpbuf[last_start + j];
                }
            }
            last_start = start;
            start += local_size;
        }
        MPI_Wait(&req_recv, &status);
        MPI_Wait(&req_send, &status);
        //开始上一轮完成传输的数据的规�
        {
#pragma omp parallel for
            for (int j = 0; j < last_local_size; j++)
            {
                recvbuf[last_start + j] += tmpbuf[last_start + j];
            }
        }

        //         MPI_Isend(recvbuf,size,MPI_CHAR,target,ALLREDUCE_TAG,allreduce_comm,&req_send);
        //         MPI_Irecv(tmpbuf,size,MPI_CHAR,target,ALLREDUCE_TAG,allreduce_comm,&req_recv);
        //         MPI_Wait(&req_recv,&status);
        //         MPI_Wait(&req_send,&status);
        // #pragma omp parallel for
        //         for(int j = 0;j<num;j++)
        //         {
        //             recvbuf[j]+=tmpbuf[j];
        //         }

        step *= 2;
    }
    free(tmpbuf);
}

extern int K_nominal_tree_stepN(int procn, int k);
extern int K_nominal_tree_my_stepN(int rank, int procn, int k);
extern int K_nominal_tree_parent(int rank, int k, int step);
extern void K_nominal_tree_child_vec(int rank, int procn, int step, int *Childvec, int *childn, int k);

#define compile_barrier() __asm__ __volatile__("" \
                                               :  \
                                               :  \
                                               : "memory");
extern MPI_Comm *vec_self_tuned_two_dimentional_allreduce_row_comm;
extern MPI_Comm *vec_self_tuned_two_dimentional_allreduce_colum_comm;
extern int vec_self_tuned_two_dimentional_allreduce_comm_NUM;
int self_tuned_tow_dimentional_allreduce_DimX = 1;
void Self_tuned_two_dimentional_allreduce_inter_node(void *sendbuf, void *recvbuf, int count, MPI_Comm row_comm, MPI_Comm colum_comm)
{
    MPI_Allreduce(sendbuf, allreduce_sendbuf, count, MPI_DOUBLE, MPI_SUM, colum_comm);
    MPI_Allreduce(allreduce_sendbuf, recvbuf, count, MPI_DOUBLE, MPI_SUM, row_comm);
    // if(global_rank == 0)
    //     printf("-----------------allreduce_sendbuf[0]=%lf\n",*((double*)allreduce_sendbuf));
    int DimX, DimY;
    int x, y;
    DimX = self_tuned_tow_dimentional_allreduce_DimX;
    int remain = inter_procn % DimX;
    MPI_Comm_rank(row_comm, &x);
    MPI_Comm_rank(colum_comm, &y);

    if (remain != 0)
    {
        DimY = inter_procn / DimX + 1;
        if (x < remain)
        {
            if (y == DimY - 1)
            {
                // receiver
                MPI_Status status;
                MPI_Recv(recvbuf, count, MPI_DOUBLE, y - 1, 0, colum_comm, &status);
            }
            else if (y == DimY - 2)
            {
                // sender
                MPI_Send(recvbuf, count, MPI_DOUBLE, y + 1, 0, colum_comm);
            }
        }
    }
}
void GLEXCOLL_Allreduce_K_ary(void *sendbuf, void *recvbuf, int startShift, int count)
{
    // puts("x");
    static struct glex_rdma_req rdma_req;
    static struct glex_rdma_req *bad_rdma_req;
    static glex_event_t *event;
    glex_ret_t ret;
    //第一步将消息拷贝到allreduce缓冲�
    memcpy(allreduce_sendbuf, sendbuf, count * sizeof(double));
    // allreduce_send_recv_pair

    if (_ReduceTreeS[_TreeID].type == LEAF)
    {
        // printf("check LEAF  %d\n",_GLEXCOLL.global_rank);
        int target = _ReduceTreeS[_TreeID].parentID;
        rdma_req.rmt_ep_addr.v = _ReduceTreeS[_TreeID].parentAddr.v;
        rdma_req.local_mh.v = send_mhs[allreduce_send_recv_pair][allreduce_rank].v;
        rdma_req.local_offset = 0;
        rdma_req.len = count * sizeof(double);
        rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
        rdma_req.rmt_offset = (count * sizeof(double)) * (_ReduceTreeS[_TreeID].DownReduceRank - 1);
        // printf("_ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n",_ReduceTreeS[_TreeID].DownReduceRank - 1);
        rdma_req.type = GLEX_RDMA_TYPE_PUT;
        rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
        rdma_req.rmt_evt.cookie[1] = 0x9696969696969696ULL;
        // rdma_req.local_evt.cookie[0] = 998;
        // rdma_req.local_evt.cookie[1] = 997;
        rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
        rdma_req.flag = GLEX_FLAG_REMOTE_EVT; //| GLEX_FLAG_LOCAL_EVT;
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
        //等待从父节点发送过来的消息
        {

            while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                ;
            if (event->cookie[1] != 997)
            {
                printf("probed a new event, but cookie[1] is invalid: %#llx\n",
                       (long long)event->cookie[1]);
            }
            PJT_discard_event(_GLEXCOLL.ep);
        }

        //将消息拷贝回用户缓冲�
        memcpy(recvbuf, allreduce_recvbuf, count * sizeof(double));
        // for(int i = 0;i<count;i++)
        // {
        //     printf("%lf ",((double *)allreduce_sendbuf)[i]);
        // }
        // puts("");
        // glex_ret_t ret;
        // glex_ep_addr_t rmt_ep_addr;
        // //puts("start wait");
        // while ((ret = glex_probe_next_mp(_GLEXCOLL.ep, &rmt_ep_addr, (void **)&(r_data), &tmp_len)) == GLEX_NO_MP)
        // {
        // }
        // for (int i = 0; i < count; i++)
        //     ((double *)recvbuf)[i] = (r_data)[i];
        // _GLEXCOLL.ep_credit -= 1;
        // if (_GLEXCOLL.ep_credit <= 0)
        // {
        //     glex_discard_probed_mp(_GLEXCOLL.ep);
        //     _GLEXCOLL.ep_credit = EP_CREDIT_MAX;
        // }

        // printf("%f \n",buf[0]);//
    }
    else if (_ReduceTreeS[_TreeID].type == MID)
    {
        {
            //接收来自children发送的消息并奖它们求和
            glex_ret_t ret;
            glex_ep_addr_t rmt_ep_addr;
            for (int i = 0; i < _ReduceTreeS[_TreeID].childsN; i++)
            {

                while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                    ;
                if (event->cookie[1] != 0x9696969696969696ULL)
                {
                    printf("probed a new event, but cookie[1] is invalid: %#llx\n",
                           (long long)event->cookie[1]);
                }
                _GLEXCOLL.event_credit -= 1; //= Event_CREDIT_MAX;
                if (_GLEXCOLL.event_credit <= 0)
                {
                    glex_discard_probed_event(_GLEXCOLL.ep);
                    _GLEXCOLL.event_credit = Event_CREDIT_MAX;
                }
                // puts("check recv");
                // printf("tmpbuf[0]=%lf buf[0]=%lf\n",tmpbuf[0],buf[0]);
            }
            for (int j = 0; j < _ReduceTreeS[_TreeID].childsN; ++j)
            {
                double *bufadd = (double *)(allreduce_recvbuf + count * j * sizeof(double));
                // printf("*bufadd=%lf\n",*bufadd);
                for (int i = 0; i < count; i++)
                {
                    ((double *)allreduce_sendbuf)[i] += bufadd[i];
                }
            }
            // printf("my rank = %d, re = %lf\n",allreduce_rank, ((double *)allreduce_sendbuf)[0]);
        }
        {
            //将收来的消息发送给父节�
            int target = _ReduceTreeS[_TreeID].parentID;
            rdma_req.rmt_ep_addr.v = _ReduceTreeS[_TreeID].parentAddr.v;
            rdma_req.local_mh.v = send_mhs[allreduce_send_recv_pair][allreduce_rank].v;
            rdma_req.local_offset = 0;
            rdma_req.len = count * sizeof(double);
            rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
            rdma_req.rmt_offset = (count * sizeof(double)) * (_ReduceTreeS[_TreeID].DownReduceRank - 1);
            // printf("_ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n",_ReduceTreeS[_TreeID].DownReduceRank - 1);
            rdma_req.type = GLEX_RDMA_TYPE_PUT;
            rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
            rdma_req.rmt_evt.cookie[1] = 0x9696969696969696ULL;
            // rdma_req.local_evt.cookie[0] = 998;
            // rdma_req.local_evt.cookie[1] = 997;
            rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
            rdma_req.flag = GLEX_FLAG_REMOTE_EVT; //| GLEX_FLAG_LOCAL_EVT;
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
        {
            //接受父节点的结果并广播给孩子
            while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                ;
            if (event->cookie[1] != 997)
            {
                printf("probed a new event, but cookie[1] is invalid: %#llx\n",
                       (long long)event->cookie[1]);
            }
            PJT_discard_event(_GLEXCOLL.ep);
            //将消息广播到孩子上去
            for (int i = 0; i < _ReduceTreeS[_TreeID].childsN; i++)
            {
                int target = _ReduceTreeS[_TreeID].childIds[i];
                rdma_req.rmt_ep_addr.v = _ReduceTreeS[_TreeID].childAddrs[i].v;
                rdma_req.local_mh.v = recv_mhs[allreduce_send_recv_pair][allreduce_rank].v;
                rdma_req.local_offset = 0;
                rdma_req.len = count * sizeof(double);
                rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
                rdma_req.rmt_offset = 0;
                rdma_req.type = GLEX_RDMA_TYPE_PUT;
                rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
                rdma_req.rmt_evt.cookie[1] = 997; //广播flag 997
                // rdma_req.local_evt.cookie[0] = 998;
                // rdma_req.local_evt.cookie[1] = 997;
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                rdma_req.flag = GLEX_FLAG_REMOTE_EVT; //| GLEX_FLAG_LOCAL_EVT;
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
            memcpy(recvbuf, allreduce_recvbuf, count * sizeof(double));
        }
    }
    else
    {
        // if( _ReduceTreeS[_TreeID].type == ROOT)
        // printf("check MID ROOT %d\n",_GLEXCOLL.global_rank);
        glex_ret_t ret;
        glex_ep_addr_t rmt_ep_addr;
        for (int i = 0; i < _ReduceTreeS[_TreeID].childsN; i++)
        {
            while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                ;
            if (event->cookie[1] != 0x9696969696969696ULL)
            {
                printf("probed a new event, but cookie[1] is invalid: %#llx\n",
                       (long long)event->cookie[1]);
            }
            _GLEXCOLL.event_credit -= 1; //= Event_CREDIT_MAX;
            if (_GLEXCOLL.event_credit <= 0)
            {
                glex_discard_probed_event(_GLEXCOLL.ep);
                _GLEXCOLL.event_credit = Event_CREDIT_MAX;
            }
            // puts("check recv");
            // printf("tmpbuf[0]=%lf buf[0]=%lf\n",tmpbuf[0],buf[0]);
        }
        for (int j = 0; j < _ReduceTreeS[_TreeID].childsN; ++j)
        {
            double *bufadd = (double *)(allreduce_recvbuf + count * j * sizeof(double));
            // printf("*bufadd=%lf\n",*bufadd);
            for (int i = 0; i < count; i++)
            {
                ((double *)allreduce_sendbuf)[i] += bufadd[i];
            }
        }
        // for(int i = 0;i<count;i++)
        // {
        //     printf("%lf ",((double *)allreduce_sendbuf)[i]);
        // }
        // puts("");
        // exit(0);
        // root要将sendbuf广播给它的孩子

        //将规约结果分发给每个孩子,
        for (int i = 0; i < _ReduceTreeS[_TreeID].childsN; i++)
        {
            int target = _ReduceTreeS[_TreeID].childIds[i];
            rdma_req.rmt_ep_addr.v = _ReduceTreeS[_TreeID].childAddrs[i].v;
            rdma_req.local_mh.v = send_mhs[allreduce_send_recv_pair][allreduce_rank].v;
            rdma_req.local_offset = 0;
            rdma_req.len = count * sizeof(double);
            rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
            rdma_req.rmt_offset = 0;
            rdma_req.type = GLEX_RDMA_TYPE_PUT;
            rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
            rdma_req.rmt_evt.cookie[1] = 997; //广播flag 997
            // rdma_req.local_evt.cookie[0] = 998;
            // rdma_req.local_evt.cookie[1] = 997;
            rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
            rdma_req.flag = GLEX_FLAG_REMOTE_EVT; //| GLEX_FLAG_LOCAL_EVT;
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
        memcpy(recvbuf, allreduce_sendbuf, count * sizeof(double));
        // //第五步是将消息拷贝到recvbuf
        // for (int i = 0; i < count; i++)
        // {
        //     ((double *)recvbuf)[i] = r_data[i];
        // }
        // if (_GLEXCOLL.ep_credit <= 0)
        // {
        //     glex_discard_probed_mp(_GLEXCOLL.ep);
        //     _GLEXCOLL.ep_credit = EP_CREDIT_MAX;
        // }
    }
    // MPI_Bcast(recvbuf,count,MPI_DOUBLE,0,allreduce_comm);
    // MPI_Barrier(Comm_inter);
}

void GLEXCOLL_Allreduce_K_ary_MPI(double *sendbuf, double *recvbuf, int startshift, int count)
{
    double *middle_re = recvbuf;
    // printf("%d send %f\n", global_rank, sendbuf[count - 1]);
    memcpy(middle_re, sendbuf, count * sizeof(double));

    if (_ReduceTreeS[_TreeID].type != LEAF)
    {
        //等待和规约
        int childn = _ReduceTreeS[_TreeID].childsN;
        MPI_Request reqs[childn];
        MPI_Status status[childn];
        for (int i = 0; i < childn; i++)
        {
            int childid = _ReduceTreeS[_TreeID].childIds[i];
            MPI_Irecv(allreduce_recvbuf + (i + 1) * count * sizeof(double), count, MPI_DOUBLE, childid, 0, Comm_inter, &(reqs[i]));
            // printf("%d recv from %d\n", inter_rank, childid);
        }
        MPI_Waitall(childn, reqs, status);

        //规约数据
        for (int i = 0; i < childn; i++)
        {
            double *startp = allreduce_recvbuf + (i + 1) * count * sizeof(double);
#pragma omp simd
            for (int j = 0; j < count; j++)
            {
                middle_re[j] += startp[j];
            }
            // printf("%d recv %f from %d\n", inter_rank, startp[count - 1], _ReduceTreeS[_TreeID].childIds[i]);
        }

        // printf("%d result = %f\n", global_rank, middle_re[count - 1]);
        // fflush(stdout);
    }
    if (_ReduceTreeS[_TreeID].type != ROOT)
    {
        //把中间结果向上发送
        MPI_Status status;
        int target = _ReduceTreeS[_TreeID].parentID;
        MPI_Send(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);
        // printf("%d senddata to %d\n", inter_rank, target);
        // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
    }
    //广播过程
    //接收数据
    if (_ReduceTreeS[_TreeID].type != ROOT)
    {
        MPI_Status status;
        int target = _ReduceTreeS[_TreeID].parentID;
        // MPI_Rend(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);
        MPI_Recv(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
        // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
    }
    if (_ReduceTreeS[_TreeID].type != LEAF)
    {
        //向下广播
        int childn = _ReduceTreeS[_TreeID].childsN;
        MPI_Request reqs[childn];
        MPI_Status status[childn];
        for (int i = 0; i < childn; i++)
        {
            int childid = _ReduceTreeS[_TreeID].childIds[i];
            MPI_Isend(middle_re, count, MPI_DOUBLE, childid, 0, Comm_inter, &(reqs[i]));
        }
        MPI_Waitall(childn, reqs, status);
    }

    // printf("%d recv %f from %d\n", inter_rank, recvbuf[count - 1]);
    // MPI_Barrier(Comm_inter);
    // exit(0);
}
// PJT使用事件避免的allreduce
//  #define EVT_AVOID
//事件避免暂时不能在立即数和RDMA之间混用，否则天河2上会产生通信速率下降问题，原因不明。
void GLEXCOLL_Allreduce_K_ary_RDMA(double *sendbuf, double *recvbuf, int startshift, int count)
{
    static struct glex_rdma_req rdma_req;
    static struct glex_rdma_req *bad_rdma_req;
    static glex_event_t *event;
    int size = count * sizeof(double) + 4;
    glex_ret_t ret;
    uint64_t k_ary_RDMA_flag_reduce = 935;
    uint64_t k_ary_RDMA_flag_bcast = 936;
    static int evt_count = 0;
    evt_count++;
    // allreduce_event_count++;
    double *middle_re = recvbuf;
    // printf("%d send %f\n", global_rank, sendbuf[count - 1]);
    memcpy(middle_re, sendbuf, count * sizeof(double));
    *(volatile int *)(allreduce_recvbuf + startshift + count * sizeof(double)) = k_ary_RDMA_flag_reduce;
    __sync_synchronize();
    if (_ReduceTreeS[_TreeID].type != LEAF)
    {
        //等待和规约
        int childn = _ReduceTreeS[_TreeID].childsN;

        for (int i = 0; i < childn; i++)
        {
#ifdef EVT_AVOID
            volatile int *waitp = allreduce_recvbuf + startshift + (i + 2) * size - 4;
            while (*waitp != k_ary_RDMA_flag_reduce)
                ;
            *waitp = 0;
#else
            {
                while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                    ;
                if (event->cookie[0] != k_ary_RDMA_flag_reduce)
                {
                    puts("GLEXCOLL_Allreduce_K_ary_RDMA 接收到未知的rdma事件，请实现事件管理机制");
                }
                // printf("childn=%d check event\n", childn);
                PJT_discard_event(_GLEXCOLL.ep);
                // volatile int *waitp = allreduce_recvbuf + startshift + (i + 2) * size - 4;
                // __sync_synchronize();
                // printf("*waitp = %d\n", *waitp);
            }
#endif
        }
        __sync_synchronize();
        //规约数据
        for (int i = 0; i < childn; i++)
        {
            double *startp = allreduce_recvbuf + startshift + (i + 1) * size;
#pragma omp simd
            for (int j = 0; j < count; j++)
            {
                middle_re[j] += startp[j];
            }
            // printf("%d recv %f from %d\n", inter_rank, startp[count - 1], _ReduceTreeS[_TreeID].childIds[i]);
        }
        // printf("%d result = %f\n", inter_rank, middle_re[count - 1]);
        // fflush(stdout);
    }
    if (_ReduceTreeS[_TreeID].type != ROOT)
    {
        //把中间结果向上发送
        int target = _ReduceTreeS[_TreeID].parentID;
        int my_down_reduce_rank = _ReduceTreeS[_TreeID].DownReduceRank;
        // MPI_Send(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);

        rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
        rdma_req.local_mh.v = recv_mhs[allreduce_send_recv_pair][allreduce_inter_rank].v;
        rdma_req.local_offset = startshift;
        rdma_req.len = size;
        rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
        rdma_req.rmt_offset = startshift + my_down_reduce_rank * size;
        // printf("%d _ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n", inter_rank, _ReduceTreeS[_TreeID].DownReduceRank - 1);
        rdma_req.type = GLEX_RDMA_TYPE_PUT;
        rdma_req.rmt_evt.cookie[0] = k_ary_RDMA_flag_reduce;
        rdma_req.rmt_evt.cookie[1] = evt_count + 1;
        // rdma_req.local_evt.cookie[0] = 996;
        // rdma_req.local_evt.cookie[1] = i + 1;
        rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
        rdma_req.flag = 0;
#ifdef EVT_AVOID
        rdma_req.flag = GLEX_FLAG_FENCE;
#else
        rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
#endif
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
        // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
    }

    //接收广播数据
    if (_ReduceTreeS[_TreeID].type != ROOT)
    {
        MPI_Status status;
        int target = _ReduceTreeS[_TreeID].parentID;
        // MPI_Rend(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);
        // MPI_Recv(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter, &status);

#ifdef EVT_AVOID
        volatile int *waitp = allreduce_recvbuf + startshift + size - 4;
        while (*waitp != k_ary_RDMA_flag_bcast)
            ;
        *waitp = 0;
#else
        {
            while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                ;
            if (event->cookie[0] != k_ary_RDMA_flag_bcast)
            {
                puts("GLEXCOLL_Allreduce_K_ary_RDMA 接收到未知的rdma事件，请实现事件管理机制");
            }
            // printf("childn=%d check event\n", childn);
            PJT_discard_event(_GLEXCOLL.ep);
        }
#endif
        // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
    }
    if (_ReduceTreeS[_TreeID].type != LEAF)
    {
        volatile int *waitp = allreduce_recvbuf + startshift + size - 4;
        // __sync_synchronize();
        *waitp = k_ary_RDMA_flag_bcast;
        __sync_synchronize();
        //向下广播
        int childn = _ReduceTreeS[_TreeID].childsN;
        for (int i = 0; i < childn; i++)
        {
            int target = _ReduceTreeS[_TreeID].childIds[i];
            {
                // MPI_Send(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);

                rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                rdma_req.local_mh.v = recv_mhs[allreduce_send_recv_pair][allreduce_inter_rank].v;
                rdma_req.local_offset = startshift;
                rdma_req.len = size;
                rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
                rdma_req.rmt_offset = startshift;
                // printf("%d _ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n", inter_rank, _ReduceTreeS[_TreeID].DownReduceRank - 1);
                rdma_req.type = GLEX_RDMA_TYPE_PUT;
                rdma_req.rmt_evt.cookie[0] = k_ary_RDMA_flag_bcast;
                rdma_req.rmt_evt.cookie[1] = evt_count + 1;
                // rdma_req.local_evt.cookie[0] = 996;
                // rdma_req.local_evt.cookie[1] = i + 1;
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                rdma_req.flag = 0;
#ifdef EVT_AVOID
                if (i == 0)
                    rdma_req.flag = GLEX_FLAG_FENCE;
                else
#else
                rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
                if (i == 0)
                    rdma_req.flag |= GLEX_FLAG_FENCE;
#endif
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
                // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
            }
            // MPI_Isend(middle_re, count, MPI_DOUBLE, childid, 0, Comm_inter, &(reqs[i]));
        }
    }

    // printf("%d result = %f\n", inter_rank, middle_re[count - 1]);
    // MPI_Barrier(Comm_inter);
    // exit(0);
    //可以开始测试正确性。
}

// #define EVT_AVOID_immRDMA

void GLEXCOLL_Allreduce_k_ary_immRDMA(double *sendbuf, double *recvbuf, int startshift, int count)
{
    static struct glex_imm_rdma_req rdma_req;
    static struct glex_imm_rdma_req *bad_rdma_req;
    static glex_event_t *event;
    glex_ret_t ret;

    int size = count * sizeof(double) + 4;
    uint64_t k_ary_RDMA_flag_reduce = 935;
    uint64_t k_ary_RDMA_flag_bcast = 936;
    static int evt_count = 0;
    evt_count++;
    // allreduce_event_count++;
    double *middle_re = recvbuf;
    // printf("%d send %f\n", global_rank, sendbuf[count - 1]);
    memcpy(middle_re, sendbuf, count * sizeof(double));
    *(volatile int *)(allreduce_recvbuf + startshift + count * sizeof(double)) = k_ary_RDMA_flag_reduce;
    __sync_synchronize();
    if (_ReduceTreeS[_TreeID].type != LEAF)
    {
        //等待和规约
        int childn = _ReduceTreeS[_TreeID].childsN;

        for (int i = 0; i < childn; i++)
        {
#ifdef EVT_AVOID_immRDMA
            volatile int *waitp = allreduce_recvbuf + startshift + (i + 2) * size - 4;
            while (*waitp != k_ary_RDMA_flag_reduce)
                ;
            *waitp = 0;
#else
            {
                while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                    ;
                if (event->cookie[0] != k_ary_RDMA_flag_reduce)
                {
                    puts("GLEXCOLL_Allreduce_K_ary_RDMA 接收到未知的rdma事件，请实现事件管理机制");
                }
                // printf("childn=%d check event\n", childn);
                PJT_discard_event(_GLEXCOLL.ep);
                // volatile int *waitp = allreduce_recvbuf + startshift + (i + 2) * size - 4;
                // __sync_synchronize();
                // printf("*waitp = %d\n", *waitp);
            }
#endif
        }
        __sync_synchronize();
        //规约数据
        for (int i = 0; i < childn; i++)
        {
            double *startp = allreduce_recvbuf + startshift + (i + 1) * size;
#pragma omp simd
            for (int j = 0; j < count; j++)
            {
                middle_re[j] += startp[j];
            }
            // printf("%d recv %f from %d\n", inter_rank, startp[count - 1], _ReduceTreeS[_TreeID].childIds[i]);
        }
        // printf("%d result = %f\n", inter_rank, middle_re[count - 1]);
        // fflush(stdout);
    }
    if (_ReduceTreeS[_TreeID].type != ROOT)
    {
        //把中间结果向上发送
        {
            //把中间结果向上发送
            int target = _ReduceTreeS[_TreeID].parentID;
            int my_down_reduce_rank = _ReduceTreeS[_TreeID].DownReduceRank;
            // MPI_Send(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);

            rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
            // rdma_req.local_mh.v = recv_mhs[allreduce_send_recv_pair][allreduce_inter_rank].v;
            // rdma_req.local_offset = startshift;
            rdma_req.data = allreduce_recvbuf + startshift;
            rdma_req.len = size;
            rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
            rdma_req.rmt_offset = startshift + my_down_reduce_rank * size;
            // printf("%d _ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n", inter_rank, _ReduceTreeS[_TreeID].DownReduceRank - 1);
            // rdma_req.type = GLEX_RDMA_TYPE_PUT;
            rdma_req.rmt_evt.cookie[0] = k_ary_RDMA_flag_reduce;
            rdma_req.rmt_evt.cookie[1] = evt_count + 1;
            // rdma_req.local_evt.cookie[0] = 996;
            // rdma_req.local_evt.cookie[1] = i + 1;
            rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
            rdma_req.flag = 0;
#ifdef EVT_AVOID_immRDMA
            rdma_req.flag = GLEX_FLAG_FENCE;
#else
            rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
#endif
            rdma_req.next = NULL;
            int ret;
            while ((ret = glex_imm_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
            {
            }
            if (ret != GLEX_SUCCESS)
            {
                if (ret == GLEX_INVALID_PARAM)
                    printf("%d, _rdma() 非法参数", global_rank);
                printf("_rdma(), return: %d\n", ret);
                exit(1);
            }
            // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
        }
        // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
    }

    //接收广播数据
    if (_ReduceTreeS[_TreeID].type != ROOT)
    {
        MPI_Status status;
        int target = _ReduceTreeS[_TreeID].parentID;
        // MPI_Rend(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);
        // MPI_Recv(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter, &status);

#ifdef EVT_AVOID_immRDMA
        volatile int *waitp = allreduce_recvbuf + startshift + size - 4;
        while (*waitp != k_ary_RDMA_flag_bcast)
            ;
        *waitp = 0;
#else
        {
            while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                ;
            if (event->cookie[0] != k_ary_RDMA_flag_bcast)
            {
                puts("GLEXCOLL_Allreduce_K_ary_RDMA 接收到未知的rdma事件，请实现事件管理机制");
            }
            // printf("childn=%d check event\n", childn);
            PJT_discard_event(_GLEXCOLL.ep);
        }
#endif
        // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
    }
    if (_ReduceTreeS[_TreeID].type != LEAF)
    {
        volatile int *waitp = allreduce_recvbuf + startshift + size - 4;
        // __sync_synchronize();
        *waitp = k_ary_RDMA_flag_bcast;
        __sync_synchronize();
        //向下广播
        int childn = _ReduceTreeS[_TreeID].childsN;
        for (int i = 0; i < childn; i++)
        {
            int target = _ReduceTreeS[_TreeID].childIds[i];
            {
                // MPI_Send(middle_re, count, MPI_DOUBLE, target, 0, Comm_inter);

                rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                // rdma_req.local_mh.v = recv_mhs[allreduce_send_recv_pair][allreduce_inter_rank].v;
                // rdma_req.local_offset = startshift;
                rdma_req.data = allreduce_recvbuf + startshift;
                rdma_req.len = size;
                rdma_req.rmt_mh.v = recv_mhs[allreduce_send_recv_pair][target].v;
                rdma_req.rmt_offset = startshift;
                // printf("%d _ReduceTreeS[_TreeID].DownReduceRank - 1 =%d\n", inter_rank, _ReduceTreeS[_TreeID].DownReduceRank - 1);
                // rdma_req.type = GLEX_RDMA_TYPE_PUT;
                rdma_req.rmt_evt.cookie[0] = k_ary_RDMA_flag_bcast;
                rdma_req.rmt_evt.cookie[1] = evt_count + 1;
                // rdma_req.local_evt.cookie[0] = 996;
                // rdma_req.local_evt.cookie[1] = i + 1;
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                rdma_req.flag = 0;
#ifdef EVT_AVOID_immRDMA
                if (i == 0)
                    rdma_req.flag = GLEX_FLAG_FENCE;
                else
#else
                rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
                if (i == 0)
                    rdma_req.flag |= GLEX_FLAG_FENCE;
#endif
                    rdma_req.next = NULL;
                int ret;
                while ((ret = glex_imm_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
                {
                }
                if (ret != GLEX_SUCCESS)
                {
                    if (ret == GLEX_INVALID_PARAM)
                        printf("%d, _rdma() 非法参数", global_rank);
                    printf("_rdma(), return: %d\n", ret);
                    exit(1);
                }
                // MPI_Recv(recvbuf, count, MPI_DOUBLE, target, 0, Comm_inter, &status);
            }
            // MPI_Isend(middle_re, count, MPI_DOUBLE, childid, 0, Comm_inter, &(reqs[i]));
        }
    }

    // printf("%d result = %f\n", inter_rank, middle_re[count - 1]);
    // MPI_Barrier(Comm_inter);
    // exit(0);
    //可以开始测试正确性。
}
void Print_K_nomial_tree_information(int rank, int procn, int k, MPI_Comm comm)
{
    int my_depth = K_nominal_tree_my_stepN(allreduce_rank, allreduce_procn, 3);
    for (int j = 0; j < procn; j++)
    {
        for (int step = 0; step < my_depth; step++)
        {
            // if (rank == 0)
            //     printf("-------------my_depth$=%d-----------------------\n", my_depth);
            if (rank == j)
            {
                int parent = K_nominal_tree_parent(rank, k, step);
                printf("step=%d myrank=%d| parent=%d\n", step, rank, parent);
                int childvec[k];
                int childn;
                K_nominal_tree_child_vec(allreduce_rank, allreduce_procn, step, childvec, &childn, k);
                if (childn > 0)
                {
                    for (int x = 0; x < childn; x++)
                    {
                        printf("------childvec[%d]=%d\n", x, childvec[x]);
                    }
                }
            }
        }
        MPI_Barrier(comm);
        MPI_Barrier(comm);
        MPI_Barrier(comm);
    }
}
// K-nominal �
int allreduce_k = 4;
extern double *allreduce_k_ary_bufvec[64];
void K_nomial_double_sum(double *sendbuf, double *recvbuf, int num)
{
    // memcpy(recvbuf, sendbuf, num * sizeof(double));
    // static MPI_Request reqVec[64];
    // static MPI_Status statusvec[64];
    // static int childvec[64];
    // // int tree_depth = K_nominal_tree_stepN(allreduce_procn, 3);
    // int k = allreduce_k;
    // int my_depth = K_nominal_tree_my_stepN(allreduce_rank, allreduce_procn, k);
    // // {
    // //     printf("%d\tmy_depth = %d\n", allreduce_rank, my_depth);
    // // }
    // //Print_K_nomial_tree_information(allreduce_rank, allreduce_procn, 3, allreduce_comm);
    // //第一步进行reduce操作
    // for (int step = 0; step < my_depth; step++)
    // {
    //     {
    //         int parent = K_nominal_tree_parent(allreduce_rank, k, step);
    //         int childvec[k];
    //         int childn;
    //         K_nominal_tree_child_vec(allreduce_rank, allreduce_procn, step, childvec, &childn, k);
    //         if (childn > 0)
    //         {
    //             //parent recv
    //             for (int c = 0; c < childn; c++)
    //             {
    //                 int child = childvec[c];
    //                 MPI_Irecv(allreduce_k_ary_bufvec[c], num, MPI_DOUBLE, child, 0, allreduce_comm, &(reqVec[c]));
    //             }
    //             for (int c = 0; c < childn; c++)
    //             {
    //                 MPI_Wait(&(reqVec[c]),&(statusvec[c]));
    //                 // for (int i = 0; i < num; i++)
    //                 //     recvbuf[i] += allreduce_k_ary_bufvec[c][i];
    //             }
    //         }
    //         if (parent != allreduce_rank)
    //         {
    //             MPI_Send(recvbuf, num, MPI_DOUBLE, parent, 0, allreduce_comm);
    //         }
    //     }
    // }
    MPI_Reduce(sendbuf, recvbuf, num, MPI_DOUBLE, MPI_SUM, 0, allreduce_comm);
    MPI_Bcast(recvbuf, num, MPI_DOUBLE, 0, allreduce_comm);
}

int using_offloading_allreduce = 1;
int GLEXCOLL_Allreduce(void *sendbuf, void *recvbuf, int count,
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    if (comm == MPI_COMM_WORLD)
    {
        allreduce_rank = global_rank;
        allreduce_procn = global_procn;
        allreduce_comm = comm;

        allreduce_intra_rank = global_rank % ppn;
        allreduce_intra_procn = ppn;
        allreduce_inter_rank = global_rank / ppn;
        allreduce_inter_procn = global_procn / ppn;

        // if (count == 1)
        // {
        //     puts("check");
        //     if (datatype == MPI_DOUBLE && op == MPI_SUM)
        //     {
        //         GLEXCOLL_Iallreduce(sendbuf, recvbuf, count, MPI_DOUBLE, MPI_SUM);
        //         GLEXCOLL_Wait_Iallreduce();
        //     }
        // }
        // else
        {
            // printf("intra_procn=%d\n",intra_procn);
            switch (SoftWare_Allreduce_Algorithm)
            {
            case Small_message_node_aware_POSIX_RDMA:
            {

                int msgsz = count * 8 + 4;
                extern int allreduce_buf_length;
                int neededL = (count * 8 + 4) * allreduce_recursive_doubling_stepn;
                static int startShift = 0;
                if (startShift + neededL >= allreduce_buf_length)
                    startShift = 0;
                if (ppn > 1)
                {
                    if (msgsz <= GLEXCOLL_max_PPSDSM_Len)
                    // if (msgsz <= 0)
                    {

                        if (using_offloading_allreduce == 1)
                        {
                            // puts("check error");
                            GLEXCOLL_Iallreduce(sendbuf, recvbuf, count, MPI_DOUBLE, MPI_SUM);
                            // usleep(uinit(e));
                            // fun(calc);
                            GLEXCOLL_Wait_Iallreduce();
                            // puts("recv");
                            return;
                        }

                        //调用PPSDSM+立即数RDMA
                        // if (global_rank == 0)
                        // {
                        //     puts("start intra reduce");
                        // }
                        // MPI_Barrier(Comm_intra);
                        one_double_allreduce_intra_REDUCE_to_0(sendbuf, allreduce_sendbuf + startShift, count);
                        // MPI_Reduce(sendbuf, allreduce_sendbuf + startShift, count, MPI_DOUBLE, MPI_SUM, 0, Comm_intra);
                        // __sync_synchronize();
                        // if (intra_rank == 0)
                        // {
                        //     printf("%d check intra reduce result=%lf\n", inter_rank, *(allreduce_sendbuf + startShift));
                        //     fflush(stdout);
                        // }
                        // MPI_Barrier(Comm_intra);
                        // {
                        //     // printf("----------------------\n reduce re: %lf\n", (double *)allreduce_sendbuf[0]);
                        // if (*(double *)(allreduce_sendbuf + startShift) - 32.0000 > 0.0001 || *(double *)(allreduce_sendbuf + startShift) - 32.0000 < -0.0001)
                        // {
                        //     puts("check error");
                        //     exit(0);
                        // }
                        // }
                        if (allreduce_inter_procn > 1)
                        {
                            double *re;
                            // static double rev[8];
                            re = recvbuf;

                            if (intra_rank == 0)
                            {
                                // printf("reduce re = %f\n", ((double *)(allreduce_sendbuf + startShift))[0]);
                                re = recursive_doubling_double_sum_imm_rdma(startShift, count);
                                // printf("bcast re = %f\n", re[0]);
                                // recursive_doubling_double_sum_rdma(allreduce_sendbuf + startShift, re, count);
                                // MPI_Allreduce(allreduce_sendbuf + startShift, re, count, MPI_DOUBLE, MPI_SUM, Comm_inter);
                                // if (global_rank == 0)
                                //     printf("re = %lf\n", *re);
                                // printf("%d check inter allreduce\n", global_rank);
                                // {
                                //     puts("check inter allreduce");
                                // }
                                one_double_allreduce_intra_BCAST_from_0(re, recvbuf, count);
                                // memcpy(recvbuf, re, sizeof(double) * count);
                                // MPI_Bcast(re, count, MPI_DOUBLE, 0, Comm_intra);
                                // if (global_rank == 0)
                                // {
                                //     puts("check intra bcast");
                                // }
                            }
                            else
                            {
                                // MPI_Bcast(recvbuf, count, MPI_DOUBLE, 0, Comm_intra);
                                one_double_allreduce_intra_BCAST_from_0(allreduce_sendbuf + startShift, recvbuf, count);
                                // if (*(double *)recvbuf - global_procn > 0.0001 || *(double *)recvbuf - global_procn < -0.0001)
                                // {
                                //     printf("check error %d %f\n", intra_rank, *(double *)recvbuf);
                                //     exit(0);
                                // }
                                // printf("%d bcast re: %lf\n", intra_rank, *(double *)recvbuf);
                            }
                        }
                        else
                        {
                            one_double_allreduce_intra_BCAST_from_0(allreduce_sendbuf + startShift, recvbuf, count);
                        }
                    }
                    else
                    {
                        //中消息allreduc
                        // puts("start");
                        double *data_recv = (double *)(allreduce_sendbuf + startShift);
                        compile_barrier();
                        medium_reduce(sendbuf, data_recv, count);
                        compile_barrier();
                        // MPI_Reduce(sendbuf, data_recv, count, datatype, op, 0, Comm_intra);
                        // if (*(double *)(allreduce_sendbuf + startShift) - intra_procn > 0.0001 || *(double *)(allreduce_sendbuf + startShift) - intra_procn < -0.0001)
                        // {
                        //     printf("nodeid = %d check error  gather result = %f\n", inter_rank, *(double *)(allreduce_sendbuf + startShift));
                        //     exit(0);
                        // }
                        // MPI_Barrier(Comm_intra);
                        if (allreduce_inter_procn > 1)
                        {
                            double *re = recvbuf;
                            if (intra_rank == 0)
                            {
                                // printf("rank=%d reduce re: %lf\n", global_rank, data_recv[0]);
                                // MPI_Allreduce(allreduce_sendbuf + startShift, recvbuf, count, datatype, op, Comm_inter);
                                recursive_doubling_double_sum_rdma(data_recv, startShift, recvbuf, count);
                                // printf("rank=%d bcast re: %lf\n", global_rank, ((double *)(recvbuf))[0]);
                                // recursive_doubling_double_sum_rdma(recvbuf, recvbuf, count);
                                // re = recursive_doubling_double_sum_rdma_avoidEvent(startShift, count);
                                // memcpy(recvbuf, re, count * sizeof(double));
                            }
                            compile_barrier();
                            medium_bcast(re, recvbuf, count);
                            // MPI_Barrier(Comm_intra);
                            // GLEX_Small_message_bcast_double(recvbuf, recvbuf, count);
                            // MPI_Bcast(recvbuf, count, datatype, 0, Comm_intra);
                            compile_barrier();
                        }
                        else
                        {
                            medium_bcast(data_recv, recvbuf, count);
                            // if (intra_rank == 0)
                            //     memcpy(recvbuf, data_recv, count * sizeof(double));
                            // MPI_Bcast(recvbuf, count, datatype, 0, Comm_intra);
                        }
                    }
                    // else
                    // {
                    //     MPI_Reduce(sendbuf, allreduce_sendbuf, count, datatype, op, 0, Comm_intra);
                    //     if (intra_rank == 0)
                    //     {
                    //         MPI_Allreduce(allreduce_sendbuf, recvbuf, count, datatype, op, Comm_inter);
                    //         memcpy(recvbuf, allreduce_recvbuf, count * sizeof(double));
                    //     }
                    //     MPI_Bcast(recvbuf, count, datatype, 0, Comm_intra);
                    // }

                    startShift += (((neededL) >> 6) << 6) + 64;
                    return 1;
                }
                else
                {
                    //每个节点单进程
                    if (msgsz <= GLEXCOLL_max_PPSDSM_Len)
                    // if (msgsz <= 0)
                    {

                        if (using_offloading_allreduce == 1)
                        {
                            // puts("check error");
                            GLEXCOLL_Iallreduce(sendbuf, recvbuf, count, MPI_DOUBLE, MPI_SUM);
                            // usleep(uinit(e));
                            // fun(calc);
                            GLEXCOLL_Wait_Iallreduce();
                            // puts("recv");
                            return;
                        }
                    }
                }
            }
            case K_nomial_tree_OMP:
                if (datatype == MPI_DOUBLE)
                {
                    int msgsz = count * 8 + 4;
                    extern int allreduce_buf_length;
                    int neededL = (count * 8 + 4) * (Childn_K + 1);
                    static int startShift = 0;
                    if (startShift + neededL >= allreduce_buf_length)
                        startShift = 0;
                    double *tmpsendbuf = allreduce_sendbuf + startShift;
                    double *tmprecvbuf = allreduce_recvbuf + startShift;
                    if (msgsz <= GLEXCOLL_max_PPSDSM_Len)
                        one_double_allreduce_intra_REDUCE_to_0(sendbuf, tmpsendbuf, count);
                    else
                        medium_reduce(sendbuf, tmpsendbuf, count);
                    // MPI_Barrier(MPI_COMM_WORLD);

                    if (intra_rank == 0)
                    {
                        // printf("reduce re = %f\n", tmpsendbuf[0]);
                        // fflush(stdout);
                        GLEXCOLL_Allreduce_K_ary_MPI(tmpsendbuf, tmprecvbuf, startShift, count);
                        // GLEXCOLL_Allreduce_K_ary_RDMA(tmpsendbuf, tmprecvbuf, startShift, count);

                        // if (count * sizeof(double) > GLEXCOLL_max_PPSDSM_Len)
                        //     GLEXCOLL_Allreduce_K_ary_RDMA(tmpsendbuf, tmprecvbuf, startShift, count);
                        // else
                        //     GLEXCOLL_Allreduce_k_ary_immRDMA(tmpsendbuf, tmprecvbuf, startShift, count);
                        // printf("bcast re = %f\n", tmprecvbuf[0]);
                    }

                    if (msgsz <= GLEXCOLL_max_PPSDSM_Len)
                        one_double_allreduce_intra_BCAST_from_0(tmprecvbuf, recvbuf, count);
                    else
                        medium_bcast(tmprecvbuf, recvbuf, count);

                    startShift += (((neededL) >> 6) << 6) + 64;
                    return 1;
                }
                break;
            case Recursize_doubling_OMP:
                if (datatype == MPI_DOUBLE)
                {
                    // puts("check 680");
                    extern void one_double_allreduce_intra_REDUCE_to_0(void *sendbuf, void *recvbuf, int count);
                    extern void one_double_allreduce_intra_BCAST_from_0(void *sendbuf, void *recvbuf, int count);
                    // if (intra_rank == 0)
                    //     puts("693");
                    if (ppn > 1)
                    {
                        one_double_allreduce_intra_REDUCE_to_0(sendbuf, recvbuf, count);
                        if (intra_rank == 0)
                            recursive_doubling_double_sum_rdma(recvbuf, 0, recvbuf, count);
                        printf("recvbuf = %f\n", *(double *)recvbuf);
                        one_double_allreduce_intra_BCAST_from_0(recvbuf, recvbuf, count);
                    }
                    else
                    {
                        recursive_doubling_double_sum_rdma(sendbuf, 0, recvbuf, count);
                    }
                    // printf("recvbuf = %f\n", *(double *)recvbuf);
                    return 1;
                }
                break;
            case Self_tuned_two_dimentional:
                if (datatype == MPI_DOUBLE)
                {
                    if (intra_procn == 1)
                    {
                        // printf("self_tuned_tow_dimentional_allreduce_DimX = %d\n",self_tuned_tow_dimentional_allreduce_DimX);
                        if (self_tuned_tow_dimentional_allreduce_DimX == 1)
                            MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, comm);
                        else
                        {
                            Self_tuned_two_dimentional_allreduce_inter_node(sendbuf, recvbuf, count,
                                                                            vec_self_tuned_two_dimentional_allreduce_row_comm[self_tuned_tow_dimentional_allreduce_DimX - 2],
                                                                            vec_self_tuned_two_dimentional_allreduce_colum_comm[self_tuned_tow_dimentional_allreduce_DimX - 2]);
                        }
                    }
                }
                return 1;
                break;
            default:
                break;
            }
            //先判断是否为没节点单进程
            if (intra_procn == 1)
            {
                //判断进程数量是否�^n个�
                if ((allreduce_procn & (allreduce_procn - 1)) == 0)
                {
                    //进程数量�^n�
                    // puts("进程数量�^n�);
                    // MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
                    switch (SoftWare_Allreduce_Algorithm)
                    {
                    case Recursize_doubling_OMP:
                        /* code */
                        if (datatype == MPI_DOUBLE)
                            recursive_doubling_double_sum(sendbuf, recvbuf, count);
                        else
                            MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
                        break;
                    case Recursize_doubling_Slicing:
                        // puts("a");
                        if (datatype == MPI_DOUBLE)
                        {
                            // puts("x");
                            recursive_doubling_double_slicing_sum(sendbuf, recvbuf, count);
                        }
                        else
                            MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
                        break;
                    default:
                        MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
                        break;
                    }
                }
                else
                {
                    //进程数量不是2^n�
                    // puts("进程数量不是2^n�);
                    MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
                }
            }
            else
            {
                //每节点多进程
                MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
            }
        }
    }
    else
    {
        MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    }
}

// ret = glex_receive_mp(_GLEXCOLL.ep, -1, &rmt_ep_addr, buf, &GLEXCOLL_databuf_len);
// if (ret == GLEX_NO_MP)
// {
// 	printf("_receive_mp() (%d), return NO_MP?\n", __LINE__);
// 	exit(1);
// }
// else if (ret != GLEX_SUCCESS)
// {
// 	printf("_receive_mp() (%d), return: %d\n", __LINE__, ret);
// 	exit(1);
// }
// while ((ret = glex_probe_next_mp(_GLEXCOLL.ep, &rmt_ep_addr,(void **) &(r_data), &tmp_len)) == GLEX_NO_MP)
// {
// }
// //printf("recv %d byte\n",tmp_len);
// _GLEXCOLL.ep_credit -= 1;
// if (_GLEXCOLL.ep_credit <= 0)
// {
//     glex_discard_probed_mp(_GLEXCOLL.ep);
//     _GLEXCOLL.ep_credit = EP_CREDIT_MAX;
// }
// GLEX_Coll_req.rmt_ep_addr.v = _ReduceTreeS[_TreeID].childAddrs[i].v;
// GLEX_Coll_req.data = r_data;
// GLEX_Coll_req.len = sizeof(double) * count;
// GLEX_Coll_req.flag = 0;
// GLEX_Coll_req.next = NULL;
// while ((ret = glex_send_imm_mp(_GLEXCOLL.ep, &GLEX_Coll_req, NULL)) == GLEX_BUSY)
// {
// }
// if (ret != GLEX_SUCCESS)
// {
//     printf("_send_imm_mp() (%d), return: %d\n", __LINE__, ret);
//     exit(1);
// }

// {
//     while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
//                 ;
//             if (event->cookie[1] != 997)
//             {
//                 printf("probed a new event, but cookie[1] is invalid: %#llx\n",
//                     (long long)event->cookie[1]);
//             }
//             _GLEXCOLL.event_credit -= 1; //= Event_CREDIT_MAX;
//             if (_GLEXCOLL.event_credit <= 0)
//             {
//                 glex_discard_probed_event(_GLEXCOLL.ep);
//                 _GLEXCOLL.event_credit = Event_CREDIT_MAX;
//             }
// }

// // if(_ReduceTreeS[_TreeID].type == ROOT)
// // 	printf("buf[0] = %f\n",tmpbuf[0]);
// //MID节点第二步是向上传输消息
// GLEX_Coll_req.rmt_ep_addr.v = _ReduceTreeS[_TreeID].parentAddr.v;
// GLEX_Coll_req.data = recvbuf;
// GLEX_Coll_req.len = sizeof(double) * count;
// GLEX_Coll_req.flag = 0;
// GLEX_Coll_req.next = NULL;
// while ((ret = glex_send_imm_mp(_GLEXCOLL.ep, &GLEX_Coll_req, NULL)) == GLEX_BUSY)
// {
// }
// if (ret != GLEX_SUCCESS)
// {
//     printf("_send_imm_mp() (%d), return: %d\n", __LINE__, ret);
//     exit(1);
// }

// //第三步是等待来自parent传输的规约结�

// while ((ret = glex_probe_next_mp(_GLEXCOLL.ep, &rmt_ep_addr, (void **)&(r_data), &tmp_len)) == GLEX_NO_MP)
// {
// }
// _GLEXCOLL.ep_credit -= 1;
// //printf("buf[0] = %f\n",buf[0]);