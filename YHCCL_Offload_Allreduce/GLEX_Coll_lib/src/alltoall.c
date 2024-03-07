#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/ipc.h>
#include <sys/stat.h>
#include <errno.h>
#include "glexcoll.h"
#include "glexalltoall.h"
#include <sys/types.h>

#define PJT_ALLTOALL_DEBUG
//整型数据压缩库
#include "vint.h"
#include "vsimple.h"
#include "vp4.h"
#include "bitpack.h"
#include "eliasfano.h"
//enum ALGORITHM {BRUCK,DIRECT};
int Alltoall_algorithm = 0;
extern int leaderN;

static int cookie[2];
MPI_Win ALLTOALL_win;
char *Bruck_buffer;
int NBUF_SEND = 4;
char SendStatus[64];
char RecvStatus[64];
char *block_tmp; //intra_procn * intra_procn * 8*1024
static char *bufS, *bufR;
uint64_t BUFSR_size;
//RDMA需要使用的相关数据
static struct glex_rdma_req rdma_req;
static struct glex_rdma_req *bad_rdma_req;
static glex_event_t *event;
static glex_mem_handle_t send_mh;

extern int leaderN;
int BruckStepN = 0; //指示着Bruck算法的运行步数

extern int mmin(int a, int b);

//复现TPDS17的工作所需要的数组和函数
unsigned int get_coordinate_x_tpds(unsigned int in);
unsigned int get_coordinate_y_tpds(unsigned int in);
int Coordinate_x[64];
int Coordinate_y[64];
int Coordinate_numa_x[4];
int Coordinate_numa_y[4];

int num_of_ongoing_msg;
int Num_of_ongoing_and_SelfAdapting_ON = 0;
int Memcpy_On = 1;

struct Runtime_Data
{
    int target;
    double time;
};
struct Runtime_Data *RT_Collections;
int Runtime_Data_cmp(const void *a, const void *b)
{
    return ((struct Runtime_Data *)a)->time > ((struct Runtime_Data *)b)->time;
}

static void glexcoll_Init_RDMA()
{
    // printf("global_rank = %d, %p\n",global_rank,_GLEXCOLL.ep);
    int ret = glex_register_mem(_GLEXCOLL.ep, bufR, BUFSR_size,
                                GLEX_MEM_READ | GLEX_MEM_WRITE,
                                &(_GLEXCOLL.local_mh));
    ret = glex_register_mem(_GLEXCOLL.ep, bufS, BUFSR_size,
                            GLEX_MEM_READ | GLEX_MEM_WRITE,
                            &send_mh);
    if (ret != GLEX_SUCCESS)
    {
        printf("_register_mem(), return: %d\n", ret);
        exit(1);
    }
    _GLEXCOLL.rmt_mhs = (glex_mem_handle_t *)malloc(sizeof(_GLEXCOLL.local_mh) * inter_procn);
    MPI_Allgather((void *)&(_GLEXCOLL.local_mh), sizeof(_GLEXCOLL.local_mh), MPI_CHAR,
                  _GLEXCOLL.rmt_mhs, sizeof(_GLEXCOLL.local_mh), MPI_CHAR, Comm_inter);
    if (intra_procn == 1)
    {
        RT_Collections = (struct Runtime_Data *)malloc(sizeof(*RT_Collections) * inter_procn);
        for (int i = 1; i < global_procn; i++)
        {

            RT_Collections[i - 1].target = (global_rank + i) % global_procn;
            RT_Collections[i - 1].time = 0.0;
        }
    }
}

static volatile void *alltoall_buffer[64];
static long long shr_seg_size;
static long long shr_buf_size = (1 << 27);
static volatile void *tpds17_alltoall_buffer;
static int tpds17_shr_seg_size;
extern void *Get_sendbuf(struct GLEXCOLL_a2a_bufmh *mh, unsigned long long i);
void glexcoll_register_alltoall_buffer(void *sendbuf, void *recvbuf, struct GLEXCOLL_a2a_bufmh *bufmh)
{
    if (ppn > 1)
    {
        uint64_t shiftsend = sendbuf - alltoall_buffer[intra_rank];
        MPI_Allgather(&shiftsend, 1, MPI_LONG, (bufmh->sendvec), 1, MPI_LONG, Comm_intra);
        uint64_t shiftrecv = recvbuf - alltoall_buffer[intra_rank];
        MPI_Allgather(&shiftrecv, 1, MPI_LONG, (bufmh->recvvec), 1, MPI_LONG, Comm_intra);
    }
    else
    {
        //当ppn == 1的时候，直接注册缓冲器的rdma地址。
        puts("每节点单进程请用 glexcoll_register_alltoall_buffer_new");
        exit(0);
    }
    // {
    //     //检查内存
    //     static int x = 0;
    //     for(int s = 0;s<ppn;s++)
    //     {
    //         char * p;
    //         for(int target = 0;target<global_procn;target++)
    //         {
    //                 p = Get_sendbuf(bufmh,s)+target *(1<<15);
    //                 x+=*p;
    //         }
    //             // for(int i = 0;i<(1<<15);i++)
    //     }
    //     MPI_Barrier(Comm_intra);
    //     if(intra_rank == 0)
    //         printf("finish mem check %d\n",x);
    // }
}

void *get_is_sendflag(int i)
{
    return alltoall_buffer[i];
}
void *get_is_recvflag(int i)
{
    return alltoall_buffer[i] + shr_buf_size + 64;
}
void *get_is_senddata_buffer(int i)
{
    return alltoall_buffer[i] + 64;
}
void *get_is_recvdata_buffer(int i)
{
    return alltoall_buffer[i] + shr_buf_size + 128;
}
// void *tpds17_get_is_sendflag(int i)
// {
//     return tpds17_alltoall_buffer + i * shr_seg_size;
// }
// void *tpds17_get_is_recvflag(int i)
// {
//     return tpds17_alltoall_buffer + i * shr_seg_size + shr_buf_size + 64;
// }
void *tpds17_get_is_senddata_buffer(int i)
{
    return tpds17_alltoall_buffer + i * tpds17_shr_seg_size;
}
void *tpds17_get_is_recvdata_buffer(int i)
{
    return tpds17_alltoall_buffer + i * tpds17_shr_seg_size + shr_buf_size;
}
void glexcoll_destroy_alltoall_shared_memory_buffer()
{
    for (int i = 0; i < intra_procn; i++)
        munmap(alltoall_buffer[i], shr_seg_size);
    munmap(tpds17_alltoall_buffer, tpds17_shr_seg_size);
}
void glexcoll_init_alltoall_shared_memory_buffer()
{
    // long long shr_buf_size_64 = 2UL*((1UL<<28)+128)
    shr_seg_size = (1UL << 28) + 128UL; //2 * (shr_buf_size + 64);
    shr_buf_size = (shr_seg_size >> 1);
    tpds17_shr_seg_size = 2 * ((1 << 14));

    int tpds17_shared_memory_size = tpds17_shr_seg_size * intra_procn;
    char name[100];
    sprintf(name, "%s-%d-alltoall\0", host_name, intra_rank);
    // {
    //     puts("142");
    //     key_t key = ftok("/tmp",intra_rank);
    //     //使用SYSTEM V共享内存
    //     int shmid = shmget(key,shr_seg_size,IPC_CREAT|0666);
    //     printf("shmid=%d\n",shmid);
    //     if(shmid < 0){ printf("shmid=%d errno=%d\n",shmid,errno);perror("shmid < 0");exit(0);}
    //     int shmids[64];
    //     MPI_Allgather(&shmid,1,MPI_INT,shmids,1,MPI_INT,Comm_intra);
    //     for (int i = 0; i < intra_procn; i++)
    //     {
    //         alltoall_buffer[i] = shmat(shmids[i],NULL,0);
    //         if(i == intra_rank)
    //             memset(alltoall_buffer[i],i,shr_seg_size);
    //         if((long)alltoall_buffer[i] == - 1) {perror("shmat error"); return ;}
    //     }
    //     MPI_Barrier(Comm_intra);
    //     puts("152");
    //     volatile int *p = get_is_senddata_buffer(intra_rank);
    //     memset(p,0,shr_buf_size);
    //     MPI_Barrier(Comm_intra);
    //     puts("156");
    //     p = get_is_recvdata_buffer(intra_rank);
    //     memset(p,0,shr_buf_size);
    //     p = (volatile int *)get_is_sendflag(intra_rank);
    //     *p = 0;
    //     p = (volatile int *)get_is_recvflag(intra_rank);
    //     *p = 0;
    // }
    // exit(0);

    {
        //每个进程分配一部分内存空间
        // printf("shr_seg_size=%d start\n",shr_seg_size);
        int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
        ftruncate64(fd, shr_seg_size);
        MPI_Barrier(Comm_intra);
        close(fd);
        // printf("shr_seg_size=%d end\n",shr_seg_size);

        //接下来要将这些内存空间映射到连续的虚拟地址上去
        for (int i = 0; i < intra_procn; i++)
        {
            sprintf(name, "%s-%d-alltoall\0", host_name, i);
            fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
            alltoall_buffer[i] = mmap64(NULL, shr_seg_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (alltoall_buffer[i] == MAP_FAILED)
            {
                perror("error mmap");
                exit(0);
            }
            if (i == intra_rank)
            {
                for (int j = 0; j < shr_seg_size; ++j)
                    ((char *)(alltoall_buffer[i]))[j] = intra_rank;
            }
            close(fd);
        }
        // puts("xxxx");
        // if(intra_rank == 0)
        // {
        //     for(int j=0;j<intra_procn;j++)
        //     {
        //         printf("%d->%d\n",((char *)(alltoall_buffer[j]))[0], ((char *)(alltoall_buffer[j]))[shr_seg_size - 1]);
        //     }
        // }
        MPI_Barrier(Comm_intra);
        // printf("*p=%d\n",((char *)(alltoall_buffer[intra_rank]))[0]);
        volatile int *p = (volatile int *)get_is_sendflag(intra_rank);
        *p = 0;
        p = (volatile int *)get_is_recvflag(intra_rank);
        *p = 0;
        volatile char *pp = (volatile char *)get_is_senddata_buffer(intra_rank);
        for (int i = 0; i < shr_buf_size; i++)
        {
            pp[i] = 0;
        }
        MPI_Barrier(Comm_intra);
        // puts("yyyy");
    }
    //tpds17所需要的共享内存
    {
        sprintf(name, "%s-tpds17alltoall\0", host_name);
        if (intra_rank == 0)
        {
            //0号进程负责分配基本空间，其它所有进程之间共享该空间
            int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
            ftruncate(fd, tpds17_shared_memory_size);
            tpds17_alltoall_buffer = mmap(NULL, tpds17_shared_memory_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            close(fd);
            volatile char *p = tpds17_alltoall_buffer;
            for (int i = 0; i < tpds17_shared_memory_size; i++)
            {
                p[i] = 0;
            }
            if (tpds17_alltoall_buffer == MAP_FAILED)
            {
                printf("error map\n");
                return;
            }
        }
        MPI_Barrier(Comm_intra);
        // puts("uuuuu");
        if (intra_rank != 0)
        {
            int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
            // ftruncate(fd, shared_memory_size);
            tpds17_alltoall_buffer = mmap(NULL, tpds17_shared_memory_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (tpds17_alltoall_buffer == MAP_FAILED)
            {
                printf("error map\n");
                return;
            }
            // printf("rank %d recv %c\n",intra_rank,*(char *)alltoall_buffer);
        }
        volatile char *p = tpds17_get_is_senddata_buffer(intra_rank);
        for (int j = 0; j < ppn; j++)
        {
            *p = 0;
            p += 64;
        }
    }
    MPI_Barrier(Comm_intra);
    // puts("check yyyy");
}

//复现TPDS17的工作，需要两个函数分别取出整数的奇数位和偶数位重新组成一个整数
//Cache-oblivious MPI all-to-all communications based on Morton order
//由于arm处理器似乎没有类似_pext_u32这样的指令。故而自己写了一个
// #include <immintrin.h>
// unsigned int get_coordinate_x(unsigned int a) {
//    return _pext_u32(a, 0x55555555);
// }

// unsigned int get_coordinate_y(unsigned int a) {
//    return _pext_u32(a, 0xAAAAAAAA);
// }

unsigned int get_coordinate_x_tpds(unsigned int in)
{
    unsigned int re = 0;
    for (int i = 30; i >= 0; i -= 2)
    {
        re |= ((in & (1 << i)) >> (i / 2));
    }
    return re;
}
unsigned int get_coordinate_y_tpds(unsigned int in)
{
    unsigned int re = 0;
    for (int i = 31; i >= 1; i -= 2)
    {
        re |= ((in & (1 << i)) >> (1 + i / 2));
    }
    return re;
}
#pragma GCC push_options
#pragma GCC optimze("O2")

void intra_memory_barrier_alltoall()
{
    __sync_synchronize();
    static int id_count;
    volatile char *get_is_buffer(int i, volatile char *p)
    {
        return p + 64 * i;
    }
    //基于PPSDSM来做一个Barrier
    volatile char *barrierp = tpds17_get_is_senddata_buffer(0);
    volatile char *p = get_is_buffer(intra_rank, barrierp);
    // printf("*p=%d\n",*p);
    // while(*p !=0);
    *p = 'R';
    // __sync_lock_test_and_set(p,'R');
    // MPI_Barrier(Comm_intra);
    if (intra_rank == 0)
    {
        //barrier 等待所有进程到达
        for (int i = 0; i < intra_procn; i++)
        {
            while (*p != 'R')
                ;
            // printf("*p=%c\n",*p);
            p += 64;
        }
        p = barrierp;
        //release 释放
        for (int i = 0; i < intra_procn; i++)
        {
            *p = 'S';
            p += 64;
        }
    }
    p = get_is_buffer(intra_rank, barrierp);
    while (*p != 'S')
        ;
    // *p = 0;
    __sync_synchronize();
}
#pragma GCC pop_options

// int TPDS17_8_coreX[8]=
void TPDS17_Cache_oblivious_intra_node_UMA_alltoall_on_buffer(int elem_size)
{
    // volatile int *  my_source_tag = get_is_sendflag(intra_rank);
    // volatile int *  my_target_tag = get_is_recvflag(intra_rank);
    //表示需要发送的数据已经准备好了
    // *my_source_tag = intra_procn;
    // printf("*my_source_tag = %d\n",*my_source_tag);
    // while(!__sync_bool_compare_and_swap((int*)my_source_tag,0,intra_procn));
    // while(!__sync_bool_compare_and_swap((int*)my_target_tag,0,intra_procn));
    intra_memory_barrier_alltoall();
    {
        for (int index = 0; index < intra_procn; index++)
        {
            int source_rank = Coordinate_x[index];
            int target_rank = Coordinate_y[index];
            // printf("source_rank=%d target_rank=%d\n", source_rank, target_rank);
            volatile char *buf_source = get_is_senddata_buffer(source_rank) + elem_size * target_rank;
            volatile char *buf_target = get_is_recvdata_buffer(target_rank) + elem_size * source_rank;
            // //等待源数据准备好
            // if(*source_tag < 0 ) puts("error source_tag");
            // while (*source_tag <= 0);
            // if(*target_tag < 0 ) puts("error target_tag");
            // while (*target_tag <= 0);
            // printf("%d->%d data = %d\n", source_rank, target_rank, *(int *)buf_source);
            for (int i = 0; i < elem_size; i++)
            {
                (buf_target)[i] = (buf_source)[i];
            }
            // memcpy(buf_target,buf_source,elem_size);
            //消息发送完成之后需要将发送tag还原为0
            //此时需要原子操作的参与

            // *source_tag -= 1;
            // __sync_fetch_and_sub((int *)source_tag, 1);
            // __sync_fetch_and_sub((int *)target_tag, 1);
        }
    }
    //进程完成时每个进程必须进行Barrier。因为每个进程需要负责其它人的通信结构。
    //当执行到此处时，还无法保证本进程所需的数据都已经准备好。
    intra_memory_barrier_alltoall();
    // MPI_Barrier(Comm_intra);
    // share_memory_barrier();
    // exit(0);
    //第二步是进行inter-cache alltoall交换
}
void TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(int shift, int elem_size)
{
    intra_memory_barrier_alltoall();
    int corePerNuma = 2;
    int numaNum = intra_procn / corePerNuma;
    int numaid = intra_rank / corePerNuma;
    int intra_numa_id = intra_rank % corePerNuma;

    for (int numaShift = 0; numaShift < numaNum; numaShift++)
    {
        int source_numa = numaid;
        int target_numa = (source_numa + numaShift) % numaNum;
        //printf("source_numa=%d target_numa=%d\n",source_numa,target_numa);
        if (corePerNuma == 2) //Matrix 3000处理器
            for (int id = 0; id < corePerNuma; id++)
            {
                int source_rank = source_numa * corePerNuma + intra_numa_id;
                int target_rank = target_numa * corePerNuma + id;
                char *buf_source = shift + get_is_senddata_buffer(source_rank) + elem_size * target_rank;
                volatile char *buf_target = shift + get_is_recvdata_buffer(target_rank) + elem_size * source_rank;
                // if(global_rank == 8)
                //     printf("source = %d,target=%d shift=%d,sourcedata=%d\n",source_rank,target_rank,shift,*(int *)buf_source);
                for (int i = 0; i < elem_size; i++)
                {
                    ((volatile char *)buf_target)[i] = ((char *)buf_source)[i];
                }
            }
        else if (corePerNuma == 4) //FT-2000+处理器
        {
            int numa_start_row = source_numa * corePerNuma;
            int numa_start_col = target_numa * corePerNuma;
            for (int i = 0; i < corePerNuma; i++)
            {
                int source_rank = numa_start_row + Coordinate_numa_y[i];
                int target_rank = numa_start_col + Coordinate_numa_x[i];
                // if(global_rank == 5)
                // printf("global_rank = %d source = %d,target=%d\n",global_rank, source_rank, target_rank);
                char *buf_source = shift + get_is_senddata_buffer(source_rank) + elem_size * target_rank;
                volatile char *buf_target = shift + get_is_recvdata_buffer(target_rank) + elem_size * source_rank;
                for (int i = 0; i < elem_size; i++)
                {
                    ((volatile char *)buf_target)[i] = ((char *)buf_source)[i];
                }
            }
            // for (int id = 0; id < 2; id++)
            // {
            //     int source_rank = source_numa * corePerNuma + (intra_numa_id & 0x2) * corePerNuma;
            //     int target_rank = target_numa * corePerNuma + id + (intra_numa_id & 0x1) * corePerNuma;
            //     char *buf_source = shift + get_is_senddata_buffer(source_rank) + elem_size * target_rank;
            //     volatile char *buf_target = shift + get_is_recvdata_buffer(target_rank) + elem_size * source_rank;
            //     if (intra_rank == 1)
            //         printf("source = %d,target=%d\n", source_rank, target_rank);
            //     for (int i = 0; i < elem_size; i++)
            //     {
            //         ((volatile char *)buf_target)[i] = ((char *)buf_source)[i];
            //     }
            // }
            // for (int id = 0; id < 2; id++)
            // {
            //     int source_rank = source_numa * corePerNuma + (intra_numa_id & 0x2) * corePerNuma + 1;
            //     int target_rank = target_numa * corePerNuma + id + (intra_numa_id & 0x1) * corePerNuma;
            //     char *buf_source = shift + get_is_senddata_buffer(source_rank) + elem_size * target_rank;
            //     volatile char *buf_target = shift + get_is_recvdata_buffer(target_rank) + elem_size * source_rank;
            //     if (intra_rank == 1)
            //         printf("source = %d,target=%d\n", source_rank, target_rank);
            //     for (int i = 0; i < elem_size; i++)
            //     {
            //         ((volatile char *)buf_target)[i] = ((char *)buf_source)[i];
            //     }
            // }
            // exit(0);
        }
    }

    intra_memory_barrier_alltoall();
}

void *Get_sendbuf(struct GLEXCOLL_a2a_bufmh *mh, unsigned long long i)
{
    return (alltoall_buffer[i] + mh->sendvec[i]);
}
void *Get_recvbuf(struct GLEXCOLL_a2a_bufmh *mh, unsigned long long i)
{
    return (alltoall_buffer[i] + mh->recvvec[i]);
}
void TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(int shift, int elem_size, struct GLEXCOLL_a2a_bufmh *bufmh)
{
    intra_memory_barrier_alltoall();
    int corePerNuma = 4;
    int numaNum = intra_procn / corePerNuma;
    int numaid = intra_rank / corePerNuma;
    int intra_numa_id = intra_rank % corePerNuma;

    for (int numaShift = 0; numaShift < numaNum; numaShift++)
    {
        int source_numa = numaid;
        int target_numa = (source_numa + numaShift) % numaNum;
        // printf("source_numa=%d target_numa=%d\n",source_numa,target_numa);
        if (corePerNuma == 2) //Matrix 3000处理器
            for (int id = 0; id < corePerNuma; id++)
            {
                int source_rank = source_numa * corePerNuma + intra_numa_id;
                int target_rank = target_numa * corePerNuma + id;
                char *buf_source = shift + Get_sendbuf(bufmh, source_rank) + elem_size * target_rank;
                volatile char *buf_target = shift + Get_recvbuf(bufmh, target_rank) + elem_size * source_rank;
                // if(global_rank == 8)
                //     printf("source = %d,target=%d shift=%d,sourcedata=%d\n",source_rank,target_rank,shift,*(int *)buf_source);
                for (int i = 0; i < elem_size; i++)
                {
                    ((volatile char *)buf_target)[i] = ((char *)buf_source)[i];
                }
            }
        else if (corePerNuma == 4) //FT-2000+处理器
        {
            int numa_start_row = source_numa * corePerNuma;
            int numa_start_col = target_numa * corePerNuma;
            for (int i = 0; i < corePerNuma; i++)
            {
                int source_rank = numa_start_row + Coordinate_numa_y[i];
                int target_rank = numa_start_col + Coordinate_numa_x[i];
                // if(global_rank == 5)
                // printf("global_rank = %d source = %d,target=%d\n",global_rank, source_rank, target_rank);
                char *buf_source = shift + Get_sendbuf(bufmh, source_rank) + elem_size * target_rank;
                volatile char *buf_target = shift + Get_recvbuf(bufmh, target_rank) + elem_size * source_rank;
                for (int i = 0; i < elem_size; i++)
                {
                    ((volatile char *)buf_target)[i] = ((char *)buf_source)[i];
                }
            }
        }
    }

    intra_memory_barrier_alltoall();
}

void shared_memory_direct_alltoall(int elem_size)
{

    // puts("check start");
    intra_memory_barrier_alltoall();
    // MPI_Barrier(Comm_intra);
    // while(__sync_bool_compare_and_swap((int*)my_target_tag,0,intra_procn));
    // puts("check");
    {
        for (int shift = 0; shift < intra_procn; shift++)
        {
            int target_rank = (intra_rank + shift) % intra_procn;
            volatile char *buf_source = get_is_senddata_buffer(intra_rank) + elem_size * target_rank;
            volatile char *buf_target = get_is_recvdata_buffer(target_rank) + elem_size * intra_rank;
            // volatile int *target_tag = get_is_recvflag(target_rank);
            //  if(*target_tag < 0 ) puts("error target_tag");
            // while (*target_tag <= 0)
            //     ;
            // printf("%d->%d data = %d\n", source_rank, target_rank, *(int *)buf_source);
            for (int i = 0; i < elem_size; i++)
            {
                ((volatile char *)buf_target)[i] = ((volatile char *)buf_source)[i];
            }
            // memcpy(buf_target,buf_source,elem_size);
            //消息发送完成之后需要将发送tag还原为0
            //此时需要原子操作的参与
            //  __sync_fetch_and_sub((int *)target_tag, 1);
            // (*(volatile int*)target_tag) -= 1;
        }
        // for (int shift = 0; shift < intra_procn; shift++)
        // {
        //     int source_rank = (intra_rank + intra_procn- shift)%intra_procn;
        //     volatile char * buf_source = get_is_senddata_buffer(source_rank) + elem_size * intra_rank;
        //     char * buf_target = get_is_recvdata_buffer(intra_rank) + elem_size * source_rank;
        //     // while (*target_tag <= 0)
        //     //     ;
        //     // printf("%d->%d data = %d\n", source_rank, target_rank, *(int *)buf_source);
        //     for (int i = 0; i < elem_size; i++)
        //     {
        //         ((char *)buf_target)[i] = ((volatile char *)buf_source)[i];
        //     }
        //     // memcpy(buf_target,buf_source,elem_size);
        //     //消息发送完成之后需要将发送tag还原为0
        //     //此时需要原子操作的参与
        //     //  __sync_fetch_and_sub((int *)target_tag, 1);
        //     //(*(volatile int*)target_tag) -= 1;
        // }
    }
    // while(__sync_bool_compare_and_swap((volatile int*)my_target_tag,intra_rank+1,intra_rank));
    // while(*my_target_tag != 0);
    // MPI_Barrier(Comm_intra);
    intra_memory_barrier_alltoall();
}
//注册用于发送接收消息的缓冲区
char *AlltoallBufSends[64];
char *AlltoallBufRecvs[64];
double alltoall_compress_ratio_sum;

void glexcoll_InitAlltoall()
{
    char *pathvar = getenv("ALLTOALL_TYPE");
    if (pathvar == 0)
    {
        puts("glexcoll_InitAlltoall pathvar=0");
        exit(1);
    }
    if (strcmp(pathvar, "BRUCK") == 0)
    {
        Alltoall_algorithm = BRUCK;
    }
    else if (strcmp(pathvar, "BRUCK_RDMA") == 0)
    {
        Alltoall_algorithm = BRUCK_RDMA;
    }
    else if (strcmp(pathvar, "DIRECT") == 0)
    {
        Alltoall_algorithm = DIRECT;
    }
    else if (strcmp(pathvar, "DIRECT_NODE_AWARE_alltoall") == 0)
    {
        Alltoall_algorithm = DIRECT_NODE_AWARE;
    }
    else if (strcmp(pathvar, "DIRECT_Kleader_NODE_AWARE") == 0)
    {
        Alltoall_algorithm = DIRECT_Kleader_NODE_AWARE;
    }
    else if (strcmp(pathvar, "DIRECT_Kleader_NODE_AWARE_RDMA") == 0)
    {
        Alltoall_algorithm = DIRECT_Kleader_NODE_AWARE_RDMA;
    }
    else if (strcmp(pathvar, "DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE") == 0)
    {
        Alltoall_algorithm = DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE;
        //节点内初始化单边通信区域。
    }
    else if (strcmp(pathvar, "XOR_EXCHANGE_RDMA_ONGOING") == 0)
    {
        Alltoall_algorithm = XOR_EXCHANGE_RDMA_ONGOING;
        //节点内初始化单边通信区域。
    }
    else if (strcmp(pathvar, "MPI_LINEAR_EXCHANGE_COMPRESS") == 0)
    {
        Alltoall_algorithm = MPI_LINEAR_EXCHANGE_COMPRESS;
        //节点内初始化单边通信区域。
    }
    else if (strcmp(pathvar, "TPDS17_Cache_oblivious_intra_node") == 0)
    {
        Alltoall_algorithm = TPDS17_Cache_oblivious_intra_node;
    }

    // puts("check 130");
    //分配缓冲区
#ifdef PJT_NEW_VERSION
    if (am_i_leader())
#else
    if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif

    {
        //BUFSR_size = (inter_procn * (1<<20) * 8 * intra_procn * intra_procn);
        block_tmp = (char *)malloc(intra_procn * intra_procn * 8 * 1024);
    }
    BUFSR_size = (1 << 25); //2147483648ULL;
    bufS = (char *)malloc(BUFSR_size);
    bufR = (char *)malloc(BUFSR_size);

    {
        AlltoallBufSends[0] = bufS; //(char *)(malloc(1<<18));
        AlltoallBufRecvs[0] = bufR; //(char *)(malloc(1<<18));
        AlltoallBufSends[1] = (char *)(malloc(BUFSR_size));
        AlltoallBufRecvs[1] = (char *)(malloc(BUFSR_size));
    }

    Bruck_buffer = (char *)malloc(BUFSR_size);
    int tmp = global_procn - 1;
    BruckStepN = 1;
    while (tmp > 1)
    {
        tmp = (tmp >> 1);
        BruckStepN++;
    }
    // puts("check 156");
    //RDMA的初始化准备
    // printf("%d %d\n",leaderN,leaderN);
#ifdef PJT_NEW_VERSION
    if (inter_procn > 1 && am_i_leader())
#else
    if (inter_procn > 1 && intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif
        glexcoll_Init_RDMA();
    //初始化节点内共享内存所需要的缓冲区
    // glexcoll_init_alltoall_shared_memory_buffer();
    //初始化TPDS17论文方法所需的转置顺序。
    {
        int start_index = intra_rank * intra_procn;
        int end_index = (intra_rank + 1) * intra_procn - 1;
        for (int i = start_index; i <= end_index; i++)
        {
            //初始化转置目标的坐标
            Coordinate_x[i - start_index] = get_coordinate_x_tpds(i);
            Coordinate_y[i - start_index] = get_coordinate_y_tpds(i);
            // if(intra_rank == 2)
            // {
            //     printf("index = %d %d %d \n",i,Coordinate_x[i-start_index],Coordinate_y[i-start_index]);
            // }
        }
    }

    //初始化NUMA-Co所需要的转置顺序
    {
        int CorePnuma = 4;
        int start_index = (intra_rank % CorePnuma) * 4;
        int end_index = start_index + 3;
        for (int i = start_index; i <= end_index; i++)
        {
            //初始化转置目标的坐标
            Coordinate_numa_x[i - start_index] = get_coordinate_x_tpds(i);
            Coordinate_numa_x[i - start_index] = get_coordinate_y_tpds(i);
            // if(intra_rank == 4)
            // {
            //     printf("global_rank=%d, index = %d %d %d \n",global_rank,i,Coordinate_x[i-start_index],Coordinate_y[i-start_index]);
            // }
        }
    }
}

void MATRIX_transform(char *buf, int n, int block_size)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < i; j++)
        {
            char *p = buf + (i * n + j) * block_size;
            char *q = buf + (j * n + i) * block_size;
            memcpy(block_tmp, p, block_size);
            memcpy(p, q, block_size);
            memcpy(q, block_tmp, block_size);
        }
    }
}
int COMMUNICATION_ON = 1;
int Intra_Gather_Scatter_ON = 1;
int MATRIX_Tanfform_ON = 1;
void DIRECT_NODE_AWARE_alltoall(void *sendbuf,
                                int sendsize,
                                void *recvbuf,
                                int recvsize,
                                MPI_Comm comm)
{
    //puts("check DIRECT_NODE_AWARE_alltoall");
    //第一步将消息聚合到leader上
    MPI_Request reqs[inter_procn];
    MPI_Request reqSend[inter_procn];
    MPI_Request reqRecv[inter_procn];
    MPI_Status status[inter_procn];
    char *source, *target;
    int block_size = sendsize * ppn;
    if (Intra_Gather_Scatter_ON)
    {
        for (int shift = 0; shift < inter_procn; ++shift)
        {
            source = sendbuf + shift * block_size;
            MPI_Gather(source, block_size, MPI_CHAR, bufS + shift * block_size * ppn, block_size, MPI_CHAR, 0, Comm_intra);
        }
    }

    if (intra_rank == 0)
    {
        if (COMMUNICATION_ON)
        {
            // MPI_Alltoall(bufS, block_size * ppn, MPI_CHAR, bufR, block_size * ppn, MPI_CHAR, Comm_inter);
            for (int shift = 0; shift < inter_procn; ++shift)
            {
                int target = (inter_rank + shift) % inter_procn;
                MPI_Isend(bufS + target * block_size * ppn, block_size * ppn, MPI_CHAR,
                          target, 0, Comm_inter, &(reqSend[shift]));
                int source = (inter_procn + inter_rank - shift) % inter_procn;
                MPI_Irecv(bufR + block_size * ppn * source, block_size * ppn, MPI_CHAR,
                          source,
                          0, Comm_inter, &(reqRecv[shift]));
            }
            MPI_Waitall(inter_procn, reqSend, status);
            MPI_Waitall(inter_procn, reqRecv, status);
        }

        //将消息转置后scatter
        if (MATRIX_Tanfform_ON)
            for (int shift = 0; shift < inter_procn; ++shift)
            {
                source = bufR + shift * block_size * ppn;
                MATRIX_transform(source, ppn, sendsize);
            }
    }
    if (Intra_Gather_Scatter_ON)
    {
        for (int shift = 0; shift < inter_procn; ++shift)
        {
            source = bufR + shift * block_size * ppn;
            target = recvbuf + shift * block_size;
            MPI_Scatter(source, block_size, MPI_CHAR, target, block_size, MPI_CHAR, 0, Comm_intra);
        }
    }
}

/*
*/

void DIRECT_Kleader_NODE_AWARE_RDMA_alltoall(void *sendbuf,
                                             int sendsize,
                                             void *recvbuf,
                                             int recvsize,
                                             MPI_Comm comm,
                                             MPI_Datatype type)
{
    MPI_Barrier(comm);
    // if(global_rank == 20)
    // {
    //     printf("%d %d\n",sendsize,recvsize);
    //     char host_name[100];
    //     int namelen, n;
    //     pid_t pid = getpid();
    //     MPI_Get_processor_name(host_name, &namelen);
    //     printf("%s %d\n",host_name,pid);
    //     volatile int pjt_d = 0;
    //     while(pjt_d == 0) ;
    // }

    MPI_Request reqs[inter_procn];
    MPI_Request reqSend[inter_procn];
    MPI_Request reqRecv[inter_procn];
    MPI_Status status[inter_procn];
    char *source, *target;
    int block_size = sendsize * ppn;
    int BufScount = 0;
    //第一步先进行对角线，节点内转置
    {
        int start_shift = (inter_rank * block_size);
        if ((sendsize >> 2) > 128)
        {
            //使用DIRECT方法转置节点内第一个大块。
            DIRECT_alltoall(sendbuf + start_shift, sendsize, recvbuf + start_shift, recvsize, Comm_intra, type);
        }
        else
        {
            //使用BMA方法转置节点内第一个大块。
            BRUCK_RMA_alltoall(sendbuf + start_shift, sendsize, recvbuf + start_shift, recvsize, Comm_intra);
        }
#ifdef PJT_ALLTOALL_DEBUG
        {
            //打印输出节点内转置结果。
            //2节点xppn核心。
            //检查正确。z
            // if(global_rank == 23)
            // {
            //     for(int i = 16;i<32;i++)
            //     {
            //         printf("%d\n",((int*)(recvbuf))[i]);
            //     }
            // }
        }
#endif
    }
    // printf("rank = %d\n",global_rank);
    //第二步进行节点间通信和转置
    {
        char *p = bufS;
        for (int shift = 1; shift < inter_procn; ++shift)
        {
            int block_id = (inter_rank + shift) % inter_procn;
            int leader = (shift % leaderN) << 2;
            // printf("%d's block_id=%d leader = %d\n",global_rank,block_id,leader);
            source = sendbuf + block_id * block_size;
            // {
            //     //打印gather输入数据
            //     if (global_rank == 16)
            //     {
            //         puts("");
            //         for(int i = 0;i<ppn;i++)
            //         {
            //             printf("gather输入数据global_rank=%d %d\n",global_rank,((int *)source)[i]);
            //         }
            //     }
            // }
            //MPI_Igather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra, &(reqs[shift]));
            MPI_Gather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra);
            {
                // //打印输出节点内gather结果。
                // if (global_rank == 4)
                // {
                //     puts("");
                //     for (int i = 0; i < ppn*ppn; i++)
                //     {
                //         printf("gather结果global_rank=%d %d\n", global_rank, ((int *)p)[i]);
                //     }
                //     puts("--------------------------------------------------");
                // }
            }
            if (intra_rank == leader)
            {
                //gather完成之后直接发送出去
                int target = (inter_rank + shift) % inter_procn;
                // printf("start send %d->%d\n",inter_rank,target); //问题，2个节点时，此处只输出start send 0->1
                rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                rdma_req.local_mh.v = send_mh.v;
                rdma_req.local_offset = (p - bufS);
                rdma_req.len = block_size * ppn;
                rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
                rdma_req.rmt_offset = BufScount * ppn * block_size;
                rdma_req.type = GLEX_RDMA_TYPE_PUT;
                rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
                rdma_req.rmt_evt.cookie[1] = 0x9696969696969696ULL;
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                if (BufScount % num_of_ongoing_msg == 1)
                    rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
                else
                    rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
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

                //更新p指向下一个gather的缓冲区。
                p += block_size * ppn;
                BufScount++;
            }
        }
#ifdef PJT_NEW_VERSION
        if (am_i_leader())
#else
        if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif
        {
            //leader节点必须等待所有消息到达
            for (int c = 0; c < BufScount; ++c)
            {
                while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                    ;
                if (event->cookie[1] != 0x9696969696969696ULL)
                {
                    printf("probed a new event, but cookie[1] is invalid: %#llx\n",
                           (long long)event->cookie[1]);
                }
            }
            //leader节点负责将消息转置
            if (MATRIX_Tanfform_ON)
                for (int c = 0; c < BufScount; ++c)
                {
                    source = bufR + c * block_size * ppn;
                    MATRIX_transform(source, ppn, sendsize);
                    // if(global_rank == 16)
                    // {
                    //     for (int i = 0; i < ppn*ppn; i++)
                    //     {
                    //         printf("转置后结果global_rank=%d %d\n", global_rank, ((int *)source)[i]);
                    //     }
                    //     exit(0);
                    // }
                }
        }
    }
    MPI_Barrier(Comm_intra);
    //第三步将所有消息scatter到对应处
    if (Intra_Gather_Scatter_ON)
    {
        char *p = bufR;
        for (int shift = 1; shift < inter_procn; ++shift)
        {
            int recv_from = (inter_procn + inter_rank - shift) % inter_procn;
            int leader = (shift % leaderN) << 2;
            // printf("%d recv_from=%d\n",global_rank,recv_from);
            // if(global_rank == 4)
            // {
            //     for (int i = 0; i < ppn*ppn; i++)
            //     {
            //         printf("Scatter前结果global_rank=%d %d\n", global_rank, ((int *)p)[i]);
            //     }
            //     puts("---------------------------------------");
            // }
            target = recvbuf + recv_from * block_size;
            MPI_Scatter(p, block_size, MPI_CHAR, target, block_size, MPI_CHAR, leader, Comm_intra);
            // if(global_rank == 4)
            // {
            //             printf("leader = %d, recv_from=%d\n",leader,recv_from);
            //             for (int i = 0; i < ppn; i++)
            //             {
            //                 printf("Scatter后结果global_rank=%d %d\n", global_rank, ((int *)target)[i]);
            //             }
            //             puts("---------------------------------------");
            // }
            if (intra_rank == leader)
                p += block_size * ppn;
        }
        //MPI_Waitall(inter_procn, reqs, status);
    }
    MPI_Barrier(comm);
    // if(global_rank == 4)
    //             for (int i = 0; i < 2*ppn; i++)
    //             {
    //                 printf("接受结果global_rank=%d %d\n", global_rank, ((int *)recvbuf)[i]);
    //             }
    //将gather过程和rdma过程流水化，原本未流水版本注释掉。
    // BufScount = inter_procn;
    // if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
    // {
    //     if (COMMUNICATION_ON)
    //     {
    //         for (int c = 1; c <= BufScount; ++c)
    //         {
    //             int shift = (intra_rank >> 2) + leaderN * c;
    //             int target = (inter_rank + shift) % inter_procn;
    //             rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
    //             // //puts("check 225");
    //             rdma_req.local_mh.v = send_mh.v;
    //             rdma_req.local_offset = c * block_size * ppn;
    //             rdma_req.len = block_size * ppn;
    //             rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
    //             rdma_req.rmt_offset = c * ppn * block_size;
    //             // //puts("check 231");
    //             rdma_req.type = GLEX_RDMA_TYPE_PUT;
    //             rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
    //             rdma_req.rmt_evt.cookie[1] = 0x9696969696969696ULL;
    //             rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
    //             if (c % num_of_ongoing_msg == 0)
    //                 rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
    //             else
    //                 rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
    //             rdma_req.next = NULL;
    //             int ret;
    //             while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
    //             {
    //             }
    //             if (ret != GLEX_SUCCESS)
    //             {
    //                 if (ret == GLEX_INVALID_PARAM)
    //                     printf("%d, _rdma() 非法参数", global_rank);
    //                 printf("_rdma(), return: %d\n", ret);
    //                 exit(1);
    //             }
    //         }
    //         for (int c = 1; c <= BufScount; ++c)
    //         {
    //             // glex_ep_addr_t tmp;
    //             // int ret = glex_receive_mp(_GLEXCOLL.ep, -1, &tmp, GLEXCOLL_databuf, &GLEXCOLL_databuf_len);
    //             // printf("recv %c\n",*(char *)&(GLEXCOLL_databuf));
    //             // glex_ep_addr_t ep_addr;
    //             // glex_get_ep_addr(_GLEXCOLL.ep, &ep_addr);
    //             // printf("local ep_addr = %#llx %#llx\n", (long long)ep_addr.v,(long long)_GLEXCOLL.ep_addrs[1].v);
    //             // printf("local ep_mem_handle = %#llx\n",
    //             //        (long long)_GLEXCOLL.local_mh.v);
    //             while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
    //                 ;
    //             if (event->cookie[1] != 0x9696969696969696ULL)
    //             {
    //                 printf("probed a new event, but cookie[1] is invalid: %#llx\n",
    //                        (long long)event->cookie[1]);
    //             }
    //         }
    //     }
    //     //将消息转置后scatter
    //     if (MATRIX_Tanfform_ON)
    //         for (int c = 1; c <= BufScount; ++c)
    //         {
    //             source = bufR + c * block_size * ppn;
    //             MATRIX_transform(source, ppn, sendsize);
    //         }
    //     glex_discard_probed_event(_GLEXCOLL.ep);
    // }

    // if (Intra_Gather_Scatter_ON)
    // {
    //     char *p = bufR;
    //     for (int shift = 0; shift < inter_procn; ++shift)
    //     {
    //         int recv_from = (inter_procn + inter_rank - shift) % inter_procn;
    //         int leader = (shift % leaderN) << 2;

    //         target = recvbuf + recv_from * block_size;
    //         MPI_Scatter(p, block_size, MPI_CHAR, target, block_size, MPI_CHAR, leader, Comm_intra);
    //         if (intra_rank == leader)
    //             p += block_size * ppn;
    //     }
    //     //MPI_Waitall(inter_procn, reqs, status);
    // }
}
int slice_size;
void DIRECT_One_slice_corePNODE_AWARE_RDMA_alltoall(void *sendbuf,
                                                    int sendsize,
                                                    void *recvbuf,
                                                    int recvsize,
                                                    MPI_Comm comm)
{
    //puts("check");
    int ret;
    int ongoin_msgN = 0;
    // int total_msgN =
    //基本策略是将每条消息分解成slice_size的大小然后流水化传输。
    //Linear exchange pattern
    int size_per_rank = sendsize * global_procn;
    int slice_count = 1;
    if (sendsize > slice_size)
        slice_count = sendsize / slice_size + ((sendsize % slice_size > 0) ? 1 : 0);

    //性能感知排序
    for (int shift = 1; shift < global_procn; shift++)
    {
        int target = (global_rank + shift) % global_procn;
        int remain_size = sendsize;
        int addr_shift = target * sendsize;
        int slice_id = 0;
        while (slice_id < slice_count)
        {
            int sizeS = mmin(remain_size, slice_size);
            //发送切片并流水化
            memcpy(bufS + addr_shift, sendbuf + addr_shift, sizeS);
            {
                rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                // if(global_rank == 1) printf("rmt_ep_addr = %lld\n",rdma_req.rmt_ep_addr.v);
                // if(global_rank == 0) printf("local_ep_addr = %lld\n",_GLEXCOLL.ep_addrs[global_rank].v);
                // //puts("check 225");
                rdma_req.local_mh.v = send_mh.v;
                rdma_req.local_offset = addr_shift;
                rdma_req.len = sizeS;
                rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
                rdma_req.rmt_offset = global_rank * sendsize + slice_id * slice_size;
                // //puts("check 231");
                rdma_req.type = GLEX_RDMA_TYPE_PUT;
                rdma_req.rmt_evt.cookie[0] = 1 + global_rank;
                rdma_req.rmt_evt.cookie[1] = 1 + slice_id;
                rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
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

            addr_shift += sizeS;
            slice_id++;
            remain_size -= sizeS;
        }
    }
    //printf("check %d\n",global_rank);
    int received_count = 0;
    while (received_count < global_procn - 1)
    {
        int received_slice = 0;
        while (received_slice < slice_count)
        {
            /* code */
            while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                ;
            int recv_from = event->cookie[0] - 1;
            int slice_id = event->cookie[1] - 1;
            int shift = recv_from * sendsize + slice_id * slice_size;
            received_slice++;
            glex_discard_probed_event(_GLEXCOLL.ep);
            int sizeTmp = slice_size;
            if (slice_id == slice_count - 1 && sendsize % slice_size != 0)
                sizeTmp = sendsize % slice_size;
            memcpy(recvbuf + shift, bufR + shift, sizeTmp);
        }
        received_count++;
    }
    memcpy(recvbuf + global_rank * sendsize, sendbuf + global_rank * sendsize, sendsize);
}

void XOR_EXCHANGE_RDMA_ONGOING_alltoall(void *sendbuf,
                                        int sendsize,
                                        void *recvbuf,
                                        int recvsize,
                                        MPI_Comm comm)
{

    //puts("check");
    int ret;
    // int total_msgN =
    //基本策略是将每条消息分解成slice_size的大小然后流水化传输。
    //Linear exchange pattern
    int size_per_rank = sendsize * global_procn;
    //memcpy(bufS,sendbuf,size_per_rank);
    //性能感知排序
    for (int shift = 1; shift < global_procn; shift++)
    {
        int target = (global_rank + shift) % global_procn;
        int remain_size = sendsize;
        int addr_shift = target * sendsize;
        if (Memcpy_On)
            memcpy(bufS + addr_shift, sendbuf + addr_shift, sendsize);
        rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
        rdma_req.local_mh.v = send_mh.v;
        rdma_req.local_offset = addr_shift;

        // for (int i = 0; i < 1; i++)
        //     printf("%d send %f\n",global_rank, ((double *)(bufS + addr_shift))[i]);
        rdma_req.len = sendsize;
        rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
        rdma_req.rmt_offset = global_rank * sendsize;
        rdma_req.type = GLEX_RDMA_TYPE_PUT;
        rdma_req.rmt_evt.cookie[0] = 1 + global_rank;
        rdma_req.rmt_evt.cookie[1] = 996; //*((double*)(bufS + addr_shift));
        rdma_req.local_evt.cookie[1] = 997;
        rdma_req.flag = GLEX_FLAG_LOCAL_EVT | GLEX_FLAG_REMOTE_EVT;
        rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
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

    //printf("check %d\n",global_rank);
    int received_count = 0;
    int sended_count = 0;
    while (received_count < global_procn - 1 || sended_count < global_procn - 1)
    {
        /* code */
        while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            ;
        int recv_from = event->cookie[0] - 1;
       cookie[1] = event->cookie[1];
        int shift = recv_from * sendsize;
        glex_discard_probed_event(_GLEXCOLL.ep);
        // for (int i = 0; i < 2; i++)
        //     printf("recv %f recv_from = %d\n", ((double *)(bufR))[i],recv_from);
        if (cookie[1] == 996)
        {
            if (Memcpy_On)
                memcpy(recvbuf + shift, bufR + shift, sendsize);
            received_count++;
        }
        else if (cookie[1] == 997)
            sended_count++;
        else
        {
            puts("error");
            exit(0);
        }

        //printf("%d check recv\n",global_rank);
    }
    //memcpy(recvbuf, bufR, global_procn*sendsize);
    if (Memcpy_On)
        memcpy(recvbuf + global_rank * sendsize, sendbuf + global_rank * sendsize, sendsize);
    MPI_Barrier(MPI_COMM_WORLD);
}

void DIRECT_One_corePNODE_AWARE_RDMA_alltoall(void *sendbuf,
                                              int sendsize,
                                              void *recvbuf,
                                              int recvsize,
                                              MPI_Comm comm)
{
    //puts("check");
    int ret;
    // int total_msgN =
    //Linear exchange pattern
    int size_per_rank = sendsize * global_procn;
    //memcpy(bufS,sendbuf,size_per_rank);
    for (int shift = 1; shift < global_procn; shift++)
    {
        int target = (global_rank + shift) % global_procn;
        int remain_size = sendsize;
        int addr_shift = target * sendsize;
        if (Memcpy_On)
            memcpy(bufS + addr_shift, sendbuf + addr_shift, sendsize);
        rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
        rdma_req.local_mh.v = send_mh.v;
        rdma_req.local_offset = addr_shift;

        // for (int i = 0; i < 1; i++)
        //     printf("%d send %f\n",global_rank, ((double *)(bufS + addr_shift))[i]);
        rdma_req.len = sendsize;
        rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
        rdma_req.rmt_offset = global_rank * sendsize;
        rdma_req.type = GLEX_RDMA_TYPE_PUT;
        rdma_req.rmt_evt.cookie[0] = 1 + global_rank;
        rdma_req.rmt_evt.cookie[1] = 996; //*((double*)(bufS + addr_shift));
        rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
        if (shift % num_of_ongoing_msg == 0)
            rdma_req.flag = GLEX_FLAG_FENCE | GLEX_FLAG_REMOTE_EVT;
        rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
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

    //printf("check %d\n",global_rank);
    int received_count = 0;
    int sended_count = 0;
    if (Memcpy_On)
        memcpy(recvbuf + global_rank * sendsize, sendbuf + global_rank * sendsize, sendsize);

    while (received_count < global_procn - 1)
    {
        /* code */
        while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            ;
        int recv_from = event->cookie[0] - 1;
       cookie[1] = event->cookie[1];
        int shift = recv_from * sendsize;
        // for (int i = 0; i < 2; i++)
        //     printf("recv %f recv_from = %d\n", ((double *)(bufR))[i],recv_from);
        // printf("%d check recv\n",global_rank);
        if (cookie[1] == 996)
        {
            if (Memcpy_On)
                memcpy(recvbuf + shift, bufR + shift, sendsize);
            received_count++;
        }
        else if (cookie[1] == 997)
            sended_count++;
        else
        {
            puts("error");
            exit(0);
        }
    }
    glex_discard_probed_event(_GLEXCOLL.ep);
    //memcpy(recvbuf, bufR, global_procn*sendsize);

    // while (glex_probe_next_event(_GLEXCOLL.ep, &event) != GLEX_NO_EVENT){
    //         glex_discard_probed_event(_GLEXCOLL.ep);
    // }
    //printf("check F%d\n",global_rank);
}
int RT_Self_Adapting_on = 0;
void DIRECT_RDMA_NUM_ONGOING_Self_Adapting(void *sendbuf,
                                           int sendsize,
                                           void *recvbuf,
                                           int recvsize,
                                           MPI_Comm comm)
{
    //每节点一进程
    //控制同时发送的消息数量
    //自适应拥塞避免
    //puts("check");

    int ret;
    int onGoingMax = num_of_ongoing_msg;
    int onGoingNum = 0;
    int received_count = 0;
    int sended_count = 0;
    //for (int shift = 1; shift < global_procn; shift++)
    for (int location = 0; location < global_procn - 1; location++)
    {
        //int shift = (global_procn + RT_Collections[location].target - global_rank)%global_procn;

        int target;

        if (RT_Self_Adapting_on)
            target = RT_Collections[location].target; //(global_rank + shift) % global_procn;//
        else
            target = (global_rank + location + 1) % global_procn;
        int remain_size = sendsize;
        int addr_shift = target * sendsize;
        memcpy(bufS + addr_shift, sendbuf + addr_shift, sendsize);
        rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
        rdma_req.local_mh.v = send_mh.v;
        rdma_req.local_offset = addr_shift;
        rdma_req.len = sendsize;
        rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
        rdma_req.rmt_offset = global_rank * sendsize;
        rdma_req.type = GLEX_RDMA_TYPE_PUT;
        rdma_req.flag = GLEX_FLAG_LOCAL_EVT | GLEX_FLAG_REMOTE_EVT;

        rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
        while (onGoingNum == onGoingMax)
        {
            //当没有空余发送信用的时候，等待发送信用
            {
                /* code */
                while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                    ;
               cookie[1] = event->cookie[1];
                if (cookie[1] == 999999999) //这是RMT EVNT
                {
                    int recv_fromORTarget = event->cookie[0] - 1;
                    int shift = recv_fromORTarget * sendsize;
                    memcpy(recvbuf + shift, bufR + shift, sendsize);
                    received_count++;
                }
                else
                {
                    double startT = event->cookie[0];
                    int location = event->cookie[1] - 1;
                    //puts("check local event");
                    sended_count++;
                    onGoingNum--;
                    double endT = MPI_Wtime();
                    if (RT_Self_Adapting_on)
                        RT_Collections[location].time = endT - startT;
                }
                glex_discard_probed_event(_GLEXCOLL.ep);
                //printf("%d check recv\n",global_rank);
            }
        }
        rdma_req.rmt_evt.cookie[0] = 1 + global_rank;
        rdma_req.rmt_evt.cookie[1] = 999999999;      //是RMT
        rdma_req.local_evt.cookie[0] = MPI_Wtime();  //请求进入队列的时间
        rdma_req.local_evt.cookie[1] = 1 + location; //007是RMT
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
        onGoingNum++;
    }
    //printf("check %d\n",global_rank);
    while (received_count < (global_procn - 1) || sended_count < (global_procn - 1))
    {
        /* code */
        while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
            ;
       cookie[1] = event->cookie[1];
        if (cookie[1] == 999999999) //这是RMT EVNT
        {
            int recv_fromORTarget = event->cookie[0] - 1;
            int shift = recv_fromORTarget * sendsize;
            memcpy(recvbuf + shift, bufR + shift, sendsize);
            received_count++;
        }
        else
        {
            double startT = event->cookie[0];
            int target = event->cookie[1] - 1;
            //puts("check local event");
            sended_count++;
            onGoingNum--;
            double endT = MPI_Wtime();
            RT_Collections[target].time = endT - startT;
        }
        glex_discard_probed_event(_GLEXCOLL.ep);
    }

    if (RT_Self_Adapting_on)
        qsort(RT_Collections, global_procn - 1, sizeof(*RT_Collections), Runtime_Data_cmp);
    memcpy(recvbuf + global_rank * sendsize, sendbuf + global_rank * sendsize, sendsize);
    MPI_Barrier(MPI_COMM_WORLD);
}

//2021/3/9
//将PIPELINE_alltoall改为支持comm输入的版本
//执行本程序必须保证comm通信子内进程数量整除节点数量。
//一个节点中只允许有一个comm通信子。
void DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE_alltoall(void *sendbuf,
                                                      int sendsize,
                                                      void *recvbuf,
                                                      int recvsize,
                                                      MPI_Comm comm)
{

    int alltoall_rank;
    int alltoall_procn;
    MPI_Comm_rank(comm, &alltoall_rank);
    MPI_Comm_size(comm, &alltoall_procn);
    // if(alltoall_rank == 0)
    // {
    //     puts("check 750");
    // }
    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Win win;
    MPI_Win_create(recvbuf, sendsize * global_procn, 1, MPI_INFO_NULL, Comm_intra, &win);
    MPI_Request reqs[inter_procn];
    MPI_Request reqSend[inter_procn];
    MPI_Request reqRecv[inter_procn];
    MPI_Status status[inter_procn];
    char *source, *target;
    int block_size = sendsize * ppn;
    int BufScount = 0;
    int received_count = 0;
    int sendsed_count = 0;
    int shift = 0;
    int shifts[100];
    char *p_send = bufS;
    char *p_recv = bufR;
    int rmt_offset = 0;
    MPI_Win_fence(0, win);
    // MPI_Barrier(MPI_COMM_WORLD);
    // puts("check 321");
    while (shift < inter_procn)
    {
        int BufScount = mmin(num_of_ongoing_msg, inter_procn - shift);
        int MySendCount = 0;
        char *p = p_send;
        if (Intra_Gather_Scatter_ON)
            for (int i = 0; i < BufScount; i++)
            {
                //现在进行shift个消息的阻塞式Gather
                int block_id = (inter_rank + shift) % inter_procn;
                int leader = (shift % leaderN) << 2;
                source = sendbuf + block_id * block_size;
                //MPI_Igather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra, &(reqs[shift]));
                MPI_Gather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra);

                if (intra_rank == leader)
                {

                    // printf("%d finished a gather \n", global_rank);
                    p += block_size * ppn;
                    shifts[MySendCount] = shift;
                    MySendCount++;
                }
                shift++;
            }
        if (COMMUNICATION_ON)
#ifdef PJT_NEW_VERSION
            if (am_i_leader())
#else
            if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif
            {

                for (int i = 0; i < MySendCount; i++)
                {
                    //现在进行shift个消息的通信请求发送
                    int target = (inter_rank + shifts[i]) % inter_procn;
                    rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                    // //puts("check 225");
                    rdma_req.local_mh.v = send_mh.v;
                    rdma_req.local_offset = p_send - bufS;
                    // if (global_rank == 32)
                    // {
                    //     for (int i = 0; i < ppn * ppn; i++)
                    //     {
                    //         printf("%d %f\n", inter_rank, ((double *)p_send)[i]);
                    //     }
                    // }

                    rdma_req.len = block_size * ppn;
                    rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
                    rdma_req.rmt_offset = rmt_offset * ppn * block_size;
                    // //puts("check 231");
                    rdma_req.type = GLEX_RDMA_TYPE_PUT;
                    rdma_req.rmt_evt.cookie[0] = 1 + inter_rank;
                    rdma_req.rmt_evt.cookie[1] = 1 + rmt_offset;
                    rmt_offset++;
                    rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                    if (i == MySendCount - 1)
                        rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
                    else
                        rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
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
                    p_send += ppn * block_size;
                    sendsed_count++;
                }
            }
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    // puts("check 393");
    if (COMMUNICATION_ON)
#ifdef PJT_NEW_VERSION
        if (am_i_leader())
#else
        if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif
        {
            while (received_count < sendsed_count)
            {
                /* code */

                while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                    ;
                received_count++;
                int recv_from = event->cookie[0] - 1;
                int offset = (event->cookie[1]) - 1;
                int shift = (inter_procn + inter_rank - recv_from) % inter_procn;

                char *received_buffer = bufR + offset * ppn * ppn * sendsize;
                // if (global_rank == 32)
                // {
                //     printf("offset =%d\n", offset);
                //     for (int i = 0; i < ppn * ppn; i++)
                //     {
                //         printf("%d %f\n", inter_rank, ((double *)bufR)[i]);
                //     }
                // }
                if (MATRIX_Tanfform_ON)
                    MATRIX_transform(received_buffer, ppn, sendsize);
                if (Intra_Gather_Scatter_ON)
                    for (int child = 0; child < ppn; child++)
                    {
                        int leader = (shift % leaderN) << 2;
                        MPI_Put(received_buffer + child * block_size, block_size, MPI_CHAR,
                                child, recv_from * ppn * sendsize, block_size, MPI_CHAR, win);
                    }
            }
            glex_discard_probed_event(_GLEXCOLL.ep);
        }

    // MPI_Barrier(MPI_COMM_WORLD);
    // puts("check 427");
    MPI_Win_fence(0, win);
    MPI_Win_free(&win);
}

//通信压缩
void MPI_LINEAR_EXCHANGE_COMPRESS_alltoall(void *sendbuf,
                                           int sendsize,
                                           void *recvbuf,
                                           int recvsize,
                                           MPI_Comm comm)
{
    MPI_Request reqs[inter_procn];
    MPI_Request reqSend[inter_procn];
    MPI_Request reqRecv[inter_procn];
    MPI_Status status[inter_procn];
    char *source, *target;
    int block_size = sendsize * ppn;
    int BufScount = 0;
    if (Intra_Gather_Scatter_ON)
    {
        char *p = bufS;
        for (int shift = 0; shift < inter_procn; ++shift)
        {
            int block_id = (inter_rank + shift) % inter_procn;
            int leader = (shift % leaderN) << 2;
            source = sendbuf + block_id * block_size;
            //MPI_Igather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra, &(reqs[shift]));
            MPI_Gather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra);

            if (intra_rank == leader)
            {
                p += block_size * ppn;
                BufScount++;
            }
        }
    }
    BufScount = inter_procn;
#ifdef PJT_NEW_VERSION
    if (am_i_leader())
#else
    if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif
    {
        if (COMMUNICATION_ON)
        {
            for (int c = 0; c < BufScount; ++c)
            {
                int shift = (intra_rank >> 2) + leaderN * c;
                int target = (inter_rank + shift) % inter_procn;
                rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
                // //puts("check 225");
                rdma_req.local_mh.v = send_mh.v;
                rdma_req.local_offset = c * block_size * ppn;
                rdma_req.len = block_size * ppn;
                rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
                rdma_req.rmt_offset = c * ppn * block_size;
                // //puts("check 231");
                rdma_req.type = GLEX_RDMA_TYPE_PUT;
                rdma_req.rmt_evt.cookie[0] = 0x9696969696969696ULL;
                rdma_req.rmt_evt.cookie[1] = 0x9696969696969696ULL;
                rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
                if (c % num_of_ongoing_msg == 0)
                    rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
                else
                    rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
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
            for (int c = 0; c < BufScount; ++c)
            {
                // glex_ep_addr_t tmp;
                // int ret = glex_receive_mp(_GLEXCOLL.ep, -1, &tmp, GLEXCOLL_databuf, &GLEXCOLL_databuf_len);
                // printf("recv %c\n",*(char *)&(GLEXCOLL_databuf));
                // glex_ep_addr_t ep_addr;
                // glex_get_ep_addr(_GLEXCOLL.ep, &ep_addr);
                // printf("local ep_addr = %#llx %#llx\n", (long long)ep_addr.v,(long long)_GLEXCOLL.ep_addrs[1].v);
                // printf("local ep_mem_handle = %#llx\n",
                //        (long long)_GLEXCOLL.local_mh.v);
                while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                    ;
                if (event->cookie[1] != 0x9696969696969696ULL)
                {
                    printf("probed a new event, but cookie[1] is invalid: %#llx\n",
                           (long long)event->cookie[1]);
                }
                glex_discard_probed_event(_GLEXCOLL.ep);
            }
        }
        //将消息转置后scatter
        if (MATRIX_Tanfform_ON)
            for (int c = 0; c < BufScount; ++c)
            {
                source = bufR + c * block_size * ppn;
                MATRIX_transform(source, ppn, sendsize);
            }
    }

    if (Intra_Gather_Scatter_ON)
    {
        char *p = bufR;
        for (int shift = 0; shift < inter_procn; ++shift)
        {
            int recv_from = (inter_procn + inter_rank - shift) % inter_procn;
            int leader = (shift % leaderN) << 2;

            target = recvbuf + recv_from * block_size;
            MPI_Scatter(p, block_size, MPI_CHAR, target, block_size, MPI_CHAR, leader, Comm_intra);
            if (intra_rank == leader)
                p += block_size * ppn;
        }
    }
}
void DIRECT_Kleader_NODE_AWARE_alltoall(void *sendbuf,
                                        int sendsize,
                                        void *recvbuf,
                                        int recvsize,
                                        MPI_Comm comm)
{
    MPI_Request reqs[inter_procn];
    MPI_Request reqSend[inter_procn];
    MPI_Request reqRecv[inter_procn];
    MPI_Status status[inter_procn];
    char *source, *target;
    int block_size = sendsize * ppn;
    int BufScount = 0;
    if (Intra_Gather_Scatter_ON)
    {
        char *p = bufS;
        for (int shift = 0; shift < inter_procn; ++shift)
        {
            int block_id = (inter_rank + shift) % inter_procn;
            int leader = (shift % leaderN) << 2;
            source = sendbuf + block_id * block_size;
            //MPI_Igather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra, &(reqs[shift]));
            MPI_Gather(source, block_size, MPI_CHAR, p, block_size, MPI_CHAR, leader, Comm_intra);

            if (intra_rank == leader)
            {
                p += block_size * ppn;
                BufScount++;
            }
        }
        // if(inter_rank == 1 && intra_rank == 4)
        // {
        //     for(int i = 0;i<ppn*ppn*BufScount;i++)
        //         printf("%f\n",((double *)bufS)[i]);
        // }
    }

#ifdef PJT_NEW_VERSION
    if (am_i_leader())
#else
    if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif
    {
        if (COMMUNICATION_ON)
        {
            //MPI_Alltoall(bufS, block_size * ppn, MPI_CHAR, bufR, block_size * ppn, MPI_CHAR, Comm_inter);
            for (int c = 0; c < BufScount; ++c)
            {
                int shift = (intra_rank >> 2) + leaderN * c;
                int target = (inter_rank + shift) % inter_procn;
                MPI_Isend(bufS + c * block_size * ppn, block_size * ppn, MPI_CHAR,
                          target, 0, Comm_inter, &(reqSend[c]));
                int source = (inter_procn + inter_rank - shift) % inter_procn;
                MPI_Irecv(bufR + c * ppn * block_size, block_size * ppn, MPI_CHAR,
                          source, 0, Comm_inter, &(reqRecv[c]));
            }
            MPI_Waitall(BufScount, reqSend, status);
            MPI_Waitall(BufScount, reqRecv, status);
        }
        // if(global_rank == 68)
        // {
        //     for(int i = 0;i<ppn*ppn*BufScount;i++)
        //         printf("%f\n",((double *)bufR)[i]);
        // }
        //将消息转置后scatter
        if (MATRIX_Tanfform_ON)
            for (int c = 0; c < BufScount; ++c)
            {
                source = bufR + c * block_size * ppn;
                MATRIX_transform(source, ppn, sendsize);
            }
        // if(global_rank == 32)
        // {
        //     for(int i = 0;i<ppn*ppn*inter_procn;i++)
        //         printf("%f\n",((double *)bufR)[i]);
        // }
        //printf("check %d\n",global_rank);
    }

    if (Intra_Gather_Scatter_ON)
    {
        char *p = bufR;
        for (int shift = 0; shift < inter_procn; ++shift)
        {
            int recv_from = (inter_procn + inter_rank - shift) % inter_procn;
            int leader = (shift % leaderN) << 2;

            target = recvbuf + recv_from * block_size;
            MPI_Scatter(p, block_size, MPI_CHAR, target, block_size, MPI_CHAR, leader, Comm_intra);
            if (intra_rank == leader)
                p += block_size * ppn;
        }
        //MPI_Waitall(inter_procn, reqs, status);
    }
}

static MPI_Request reqSV[32];
static MPI_Request reqRV[32];
static MPI_Status statusV[32];
//适用于每节点任意进程
//通信压缩启用
void DIRECT_Compress_alltoall(void *sendbuf,
                              int sendsize,
                              void *recvbuf,
                              int recvsize,
                              MPI_Comm comm,
                              MPI_Datatype type)
{
    // MPI_Request req_send, req_recv;
    // MPI_Status status;
    // int group_rank;
    // int group_procn;
    // MPI_Comm_rank(comm, &group_rank);
    // MPI_Comm_size(comm, &group_procn);
    // memcpy(recvbuf + sendsize * group_rank,
    //        sendbuf + recvsize * group_rank,
    //        sendsize);
    // int source, target;
    // int prev_S, prev_T;
    // int shift = 1;
    // int bufid = shift & 1;
    // int prevbufid = bufid;
    // alltoall_compress_ratio_sum = 2.0;
    // {
    //     //对于shift=1时，把消息直接发送出去，不做任何压缩
    //     {
    //         target = (group_rank + shift) % group_procn;
    //         source = (group_rank + group_procn - shift) % group_procn;
    //         prev_S = source;
    //         prev_T = target;

    //         MPI_Isend(sendbuf + sendsize * target, sendsize, MPI_CHAR,
    //                   target, 0, comm, &req_send);
    //         MPI_Irecv(AlltoallBufRecvs[bufid], recvsize, MPI_CHAR, source, 0, comm, &req_recv);
    //     }
    //     // puts("1122");
    //     for (shift = 2; shift < group_procn; shift++)
    //     {
    //         // puts("1124");
    //         int compressed_size = 0;
    //         {
    //             //压缩当前步消息
    //             prevbufid = bufid;
    //             bufid = shift & 1;
    //             prev_S = source;
    //             prev_T = target;
    //             target = (group_rank + shift) % group_procn;
    //             source = (group_rank + group_procn - shift) % group_procn;
    //             if (type == MPI_INT)
    //             {
    //                 compressed_size = p4nenc128v32((int *)(sendbuf + sendsize * target),
    //                                                (sendsize >> 2),
    //                                                (int *)AlltoallBufSends[bufid]);
    //                 // if(group_rank == 0)
    //                 // {
    //                 //     printf("%lf\n",compressed_size/sendsize);
    //                 // }
    //                 if (compressed_size > sendsize)
    //                 {
    //                     puts("check error compressed_size > sendsize");
    //                     exit(0);
    //                 }
    //                 alltoall_compress_ratio_sum += (compressed_size / sendsize);
    //             }
    //         }
    //         {
    //             MPI_Wait(&req_send, statusV);
    //             MPI_Wait(&req_recv, statusV);
    //             // puts("check start");
    //             MPI_Isend(AlltoallBufSends[bufid], compressed_size, MPI_CHAR, target, 0, comm, &req_send);
    //             MPI_Irecv(AlltoallBufRecvs[bufid], recvsize, MPI_CHAR, source, 0, comm, &req_recv);
    //         }
    //         if (shift != 2)
    //         {
    //             //将上一步消息解压缩到 recvbuf + alltoall_single_data_size*rdispls[source]
    //             // p4ndec128v32((int *)AlltoallBufRecvs[prevbufid],
    //             //              (recvsize >> 2),
    //             //              recvbuf + recvsize * prev_S);
    //         }
    //         else
    //         {
    //             memcpy(recvbuf + recvsize * prev_S,
    //                    AlltoallBufRecvs[prevbufid],
    //                    recvsize);
    //         }
    //     }

    //     {
    //         MPI_Wait(&req_send, statusV);
    //         MPI_Wait(&req_recv, statusV);
    //         // puts("check");
    //         // if (group_procn > 2)
    //         //     p4ndec128v32((int *)AlltoallBufRecvs[bufid],
    //         //                  (recvsize >> 2),
    //         //                  recvbuf + recvsize * source);
    //         // else
    //         // {
    //         //     memcpy(recvbuf + recvsize * source,
    //         //            AlltoallBufRecvs[bufid],
    //         //            recvsize);
    //         // }
    //     }
    // }
    // if(group_rank == 10)
    // printf("alltoall_compress_ratio_sum = %lf\n",alltoall_compress_ratio_sum/group_procn);
}
void DIRECT_alltoall(void *sendbuf,
                     int sendsize,
                     void *recvbuf,
                     int recvsize,
                     MPI_Comm comm,
                     MPI_Datatype type)
{
    //MPI_Alltoall(sendbuf,sendsize,MPI_DOUBLE,recvbuf,recvsize,MPI_DOUBLE,comm);
    // MPI_Request req_sends[global_procn];
    // MPI_Request req_recvs[global_procn];
    // MPI_Status status[global_procn];
    int alltoall_rank;
    int alltoall_procn;
    MPI_Comm_rank(comm, &alltoall_rank);
    MPI_Comm_size(comm, &alltoall_procn);
    // if(global_rank ==1)
    //     printf("sendbuf[0] = %f,sendbuf[1] = %f\n",((double *)sendbuf)[0],((double *)sendbuf)[1]);
    for (int shift_start = 0; shift_start < alltoall_procn; shift_start += num_of_ongoing_msg)
    {
        int s;
        for (s = 0; s < num_of_ongoing_msg && (s + shift_start) < alltoall_procn; s++)
        {
            int shift = shift_start + s;
            {
                int target = (alltoall_rank + shift) % alltoall_procn;
                MPI_Isend(sendbuf + sendsize * target, sendsize, MPI_CHAR,
                          target, alltoall_rank, comm, &(reqSV[s]));
                int source = (alltoall_procn + alltoall_rank - shift) % alltoall_procn;
                MPI_Irecv(recvbuf + recvsize * source, sendsize, MPI_CHAR,
                          source,
                          source, comm, &(reqRV[s]));
                // MPI_Wait(&(req_sends[shift]), status);
                // MPI_Wait(&(req_recvs[shift]), status);
            }
        }
        MPI_Waitall(s, reqSV, statusV);
        MPI_Waitall(s, reqRV, statusV);
    }
}

void BRUCK_RMA_alltoall(void *sendbuf,
                        int sendsize,
                        void *recvbuf,
                        int recvsize,
                        MPI_Comm comm)
{
    int alltoall_rank;
    int alltoall_procn;
    MPI_Comm_rank(comm, &alltoall_rank);
    MPI_Comm_size(comm, &alltoall_procn);

    int data_send = 0;
    // char * start = sendbuf + global_rank*sendsize;
    char *tmpbuf = (char *)malloc(sizeof(char) * alltoall_procn * sendsize);
    for (int i = 0; i < alltoall_procn; i++)
    {
        char *source = sendbuf + sendsize * ((alltoall_rank + i) % alltoall_procn);
        char *target = tmpbuf + sendsize * i;
        memcpy(target, source, sendsize);
    }
    int tmp = alltoall_procn;
    char *bufPack = (char *)malloc(sizeof(char) * alltoall_procn * sendsize);
    // if(global_rank == 1)
    //     printf("BruckStepN= %d\n",BruckStepN);
    // if(global_rank==2){
    //     printf("BEFORE: tmpbuf[0]= %f tmpbuf[1]=%f tmpbuf[2]=%f\n",((double*)tmpbuf)[0],((double*)tmpbuf)[1],((double*)tmpbuf)[2]);
    // }
    // if( global_rank == 0)
    // {
    //     printf("global rank = %d, BruckStepN =%d\n",global_rank,BruckStepN);
    // }
    for (int step = 0; step < BruckStepN; step++)
    {
        char *bufPackP = bufPack;
        unsigned int checkValue = (1 << step);
        // if (global_rank == 1)
        //     printf("-----checkValue = %d -------\n",checkValue);
        //接下来将消息打包到bufPack
        for (unsigned int shift = 1; shift < alltoall_procn; ++shift)
        {
            // if (global_rank == 1)
            //     printf("shift = %d\n",shift);
            if ((shift & checkValue))
            {
                char *source = tmpbuf + sendsize * shift;
                memcpy(bufPackP, source, sendsize);
                bufPackP += sendsize;
                // if (global_rank == 1)
                // {
                //     printf("pack %f  ", ((double *)bufPack)[0]);
                // }
            }
        }
        if (COMMUNICATION_ON)
        {
            MPI_Request req;
            MPI_Status status;
            //消息打包完毕后,开始数据传输。
            // MPI_Put(bufPack, (bufPackP - bufPack), MPI_CHAR,
            //         (global_rank + checkValue) % global_procn,
            //         0,
            //         (bufPackP - bufPack),
            //         MPI_CHAR,
            //         ALLTOALL_win);
            MPI_Irecv(Bruck_buffer, (bufPackP - bufPack), MPI_CHAR, (alltoall_procn + alltoall_rank - checkValue) % alltoall_procn, 0, comm, &req);
            MPI_Send(bufPack, (bufPackP - bufPack), MPI_CHAR, (alltoall_rank + checkValue) % alltoall_procn, 0, comm);
            MPI_Wait(&req, &status);
        }
        //每一步传输完成后开始消息解包
        char *bufUnpackP = Bruck_buffer;
        for (int shift = 1; shift < alltoall_procn; ++shift)
        {
            if ((shift & checkValue))
            {
                char *target = tmpbuf + sendsize * shift;
                memcpy(target, bufUnpackP, recvsize);
                bufUnpackP += recvsize;
            }
        }
    }
    // if (global_rank == 2)
    // {
    //     printf("tmpbuf[0]= %f tmpbuf[1]=%f tmpbuf[2]=%f\n", ((double *)tmpbuf)[0], ((double *)tmpbuf)[1], ((double *)tmpbuf)[2]);
    // }
    //最后将数据还原到recvbuf中去
    for (int i = 0; i < alltoall_procn; i++)
    {
        char *target = recvbuf + sendsize * i;
        char *source = tmpbuf + sendsize * ((alltoall_procn + alltoall_rank - i) % alltoall_procn);
        memcpy(target, source, sendsize);
    }
    // if (global_rank == 2)
    // {
    //     printf("recvbuf[0]= %f recvbuf[1]=%f recvbuf[2]=%f\n", ((double *)recvbuf)[0], ((double *)recvbuf)[1], ((double *)recvbuf)[2]);
    // }

    free(tmpbuf);
    free(bufPack);
}

void BRUCK_RDMA_alltoall(void *sendbuf,
                         int sendsize,
                         void *recvbuf,
                         int recvsize,
                         MPI_Comm comm)
{
    int data_send = 0;
    char *tmpbuf = (char *)malloc(sizeof(char) * global_procn * sendsize);
    for (int i = 0; i < global_procn; i++)
    {
        char *source = sendbuf + sendsize * ((global_rank + i) % global_procn);
        char *target = tmpbuf + sendsize * i;
        if (Memcpy_On)
            memcpy(target, source, sendsize);
    }
    int tmp = global_procn;
    int recv_flag[20];
    memset(recv_flag, 0, sizeof(int) * 20);
    for (int step = 0; step < BruckStepN; step++)
    {
        int memshift = sendsize * global_procn * step;
        char *bufPackP = bufS + memshift;
        unsigned int checkValue = (1 << step);
        for (unsigned int shift = checkValue; shift < global_procn; shift += (checkValue << 1))
        {
            int lenb = sendsize * mmin(checkValue, global_procn - shift);
            if (Memcpy_On)
                memcpy(bufPackP, tmpbuf + sendsize * shift, lenb);
            bufPackP += lenb;
        }
        if (COMMUNICATION_ON)
        {
            // MPI_Request req;
            // MPI_Status status;
            // MPI_Irecv(bufR, (bufPackP - bufS), MPI_CHAR, (global_procn + global_rank - checkValue) % global_procn, 0, MPI_COMM_WORLD, &req);
            // MPI_Send(bufS,(bufPackP - bufS), MPI_CHAR, (global_rank + checkValue) % global_procn, 0, MPI_COMM_WORLD);
            // MPI_Wait(&req, &status);

            int target = (global_rank + checkValue) % global_procn;
            rdma_req.rmt_ep_addr.v = _GLEXCOLL.ep_addrs[target].v;
            rdma_req.local_mh.v = send_mh.v;
            rdma_req.local_offset = memshift;
            // for (int i = 0; i < 1; i++)
            //     printf("%d send %f\n",global_rank, ((double *)(bufS + addr_shift))[i]);
            rdma_req.len = (bufPackP - bufS - memshift);
            rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
            rdma_req.rmt_offset = memshift;
            rdma_req.type = GLEX_RDMA_TYPE_PUT;
            rdma_req.rmt_evt.cookie[0] = 1 + global_rank;
            rdma_req.rmt_evt.cookie[1] = 1 + step;   //*((double*)(bufS + addr_shift));
            rdma_req.local_evt.cookie[1] = 1 + step; //*((double*)(bufS + addr_shift));
            rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
            rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
            int ret;
            while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
                ;
            {
            }
            if (ret != GLEX_SUCCESS)
            {
                if (ret == GLEX_INVALID_PARAM)
                    printf("%d, _rdma() 非法参数", global_rank);
                printf("_rdma(), return: %d\n", ret);
                exit(1);
            }
            // int recv_from = event->cookie[0] - 1;
            // int shift = recv_from * sendsize;
            //等待传输完成
        }
        while (recv_flag[step] != 1)
        {
            /* code */
            while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
                ;
            int step_rmt = event->cookie[1] - 1;
            //printf("%d step_rmt = %d\n",global_rank,step_rmt);
            recv_flag[step_rmt] = 1;
            glex_discard_probed_event(_GLEXCOLL.ep);
        }

        //每一步传输完成后开始消息解包
        char *bufUnpackP = bufR + memshift;
        for (int shift = 1; shift < global_procn; ++shift)
        {
            if ((shift & checkValue))
            {
                char *target = tmpbuf + sendsize * shift;
                if (Memcpy_On)
                    memcpy(target, bufUnpackP, recvsize);
                bufUnpackP += recvsize;
            }
        }
    }
    //最后将数据还原到recvbuf中去

    for (int i = 0; i < global_procn; i++)
    {
        char *target = recvbuf + sendsize * i;
        char *source = tmpbuf + sendsize * ((global_procn + global_rank - i) % global_procn);
        if (Memcpy_On)
            memcpy(target, source, sendsize);
    }
    free(tmpbuf);
    MPI_Barrier(MPI_COMM_WORLD);
}

int GLEXCOLL_Alltoall(void *sendbuf,
                      int sendsize,
                      void *recvbuf,
                      int recvsize,
                      MPI_Comm comm,
                      MPI_Datatype type)
{
    switch (Alltoall_algorithm)
    {
    case DIRECT:
        /* code */
        DIRECT_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm, type);
        break;
    case BRUCK:
        BRUCK_RMA_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
        break;
    case BRUCK_RDMA:
        BRUCK_RDMA_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
        break;
    case DIRECT_NODE_AWARE:
        if (intra_procn != 1)
            DIRECT_NODE_AWARE_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
        break;
    case DIRECT_Kleader_NODE_AWARE:
        if (intra_procn != 1)
            DIRECT_Kleader_NODE_AWARE_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
        break;
    case DIRECT_Kleader_NODE_AWARE_RDMA:
        //if(global_rank == 0) printf("intra_procn = %d\n",intra_procn);
        if (intra_procn != 1)
            DIRECT_Kleader_NODE_AWARE_RDMA_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm, type);
        else if (intra_procn == 1)
        {
            if (Num_of_ongoing_and_SelfAdapting_ON)
            {
                DIRECT_RDMA_NUM_ONGOING_Self_Adapting(sendbuf, sendsize, recvbuf, recvsize, comm);
            }
            else
            {
                //puts("check");
                DIRECT_One_corePNODE_AWARE_RDMA_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
            }
        }
        //DIRECT_One_slice_corePNODE_AWARE_RDMA_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
        break;
    case DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE:
        if (inter_procn != 1)
            DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
        break;
    case MPI_LINEAR_EXCHANGE_COMPRESS:
        if ((sendsize >> 2) < 128)
        {
            DIRECT_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm, type);
        }
        else
        {
            DIRECT_Compress_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm, type);
        }
        break;

    case XOR_EXCHANGE_RDMA_ONGOING:
        if (intra_procn == 1)
            XOR_EXCHANGE_RDMA_ONGOING_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
        break;
    case TPDS17_Cache_oblivious_intra_node:
        //一定要记得加break，不然约测越慢
        TPDS17_Cache_oblivious_intra_node_UMA_alltoall_on_buffer(sendsize);
        break;
    case shared_memory_direct_intra_node:
        //一定要记得加break，不然约测越慢
        shared_memory_direct_alltoall(sendsize);
        break;
    case TPDS17_Cache_oblivious_intra_node_NUMA:
        //一定要记得加break，不然约测越慢
        TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(0, sendsize);
        break;
    default:
        break;
    }
}

void GLEXCOLL_AlltoallFinalize()
{
    free(Bruck_buffer);
    if (intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
    {
        free(bufS);
        free(bufR);
        free(block_tmp);
    }
    //
    if (intra_rank == 0)
        free(RT_Collections);

        //由于现在GLEX FAST端口不足，无法让每个进程开启一个FAST端口
#ifdef PJT_NEW_VERSION
    if (inter_procn > 1 && am_i_leader())
#else
    if (inter_procn > 1 && intra_rank % 4 == 0 && (intra_rank >> 2) < leaderN)
#endif
    {

        int ret = glex_deregister_mem(_GLEXCOLL.ep, _GLEXCOLL.local_mh);
        if (ret != GLEX_SUCCESS)
        {
            printf("_deregister_mem(), return: %d\n", ret);
            exit(1);
        }

        ret = glex_deregister_mem(_GLEXCOLL.ep, send_mh);
        if (ret != GLEX_SUCCESS)
        {
            printf("_deregister_mem(), return: %d\n", ret);
            exit(1);
        }
    }
}

// int main()
// {
//     for(unsigned int i = 0;i<64;i++)
//     {
//         printf("%d's x=%d %d y=%d %d\n",i,get_coordinate_x(i),get_coordinate_x_tpds(i),get_coordinate_y(i),get_coordinate_y_tpds(i));
//     }
//     return 0;
// }