
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <pthread.h>
#include "glexcoll.h"

int numa_noden;

int bufN_between_ab = 128;
int bufsize = (1 << 20);

struct SHM_flag_between_source_dest
{
    char flag[64];
    int size[64];
    volatile int front;
    volatile int tail;
};

static volatile struct SHM_flag_between_source_dest *SHM_flags_Recv[64];
static volatile struct SHM_flag_between_source_dest *SHM_flags_Send[64];
int Using_buf_id[64]; //Using_buf_id[source][dest] in [0,bufN_between_ab)

//每对进程间的消息大小bufN_between_ab * bufsize
void Init_SHMBUF()
{
    int size = intra_procn * bufN_between_ab * bufsize;
    char name[100];
    sprintf(name, "%s-%d\0", host_name, intra_rank);
    int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
    ftruncate(fd, size);
    volatile char *p = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (p == MAP_FAILED)
    {
        printf("error map\n");
        return;
    }
    memset((void *)p, 'S' + intra_rank, size);
    for (int i = 0; i < intra_procn; ++i)
    {

        SHM_bufs_Recv[i] = p + i * bufN_between_ab * bufsize;
        ((int *)(SHM_bufs_Recv[i]))[0] = i;
        ((int *)(SHM_bufs_Recv[i]))[1] = intra_rank;
    }
    MPI_Barrier(Comm_intra);
    //接下来记录下每对进程之间共享内存的区域
    for (int i = 0; i < intra_procn; ++i)
    {
        sprintf(name, "%s-%d\0", host_name, i);
        fd = shm_open(name, O_RDWR, S_IRUSR | S_IWUSR);
        ftruncate(fd, size);
        SHM_bufs_Send[i] = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        SHM_bufs_Send[i] += intra_rank * bufN_between_ab * bufsize;
        // if(intra_rank == 1)
        //     printf("%d -> %d \n",((int*)(SHM_bufs_Send[i]))[0],((int*)(SHM_bufs_Send[i]))[1] );
    }
    MPI_Barrier(Comm_intra);
}

void Init_SHEMFLAG()
{
    int size = intra_procn * sizeof(struct SHM_flag_between_source_dest);
    char name[100];
    sprintf(name, "%s-%d-FLAG\0", host_name, intra_rank);
    int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
    ftruncate(fd, size);
    volatile char *p = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (p == MAP_FAILED)
    {
        printf("error map\n");
        return;
    }
    for (int i = 0; i < intra_procn; ++i)
    {

        SHM_flags_Recv[i] = (volatile struct SHM_flag_between_source_dest *)(p + i * sizeof(struct SHM_flag_between_source_dest));
        for (int j = 0; j < bufN_between_ab; ++j)
        {
            SHM_flags_Recv[i]->flag[j] = 'S';
            SHM_flags_Recv[i]->size[j] = 0;
            SHM_flags_Recv[i]->front = 0;
            SHM_flags_Recv[i]->tail = 0;
        }
        //    ((int *)(SHM_flags_Recv[i]))[0] = i;
        //    ((int *)(SHM_flags_Recv[i]))[1] = intra_rank;
    }
    MPI_Barrier(Comm_intra);
    //puts("check");
    //接下来记录下每对进程之间共享内存的区域
    for (int i = 0; i < intra_procn; ++i)
    {
        sprintf(name, "%s-%d-FLAG\0", host_name, i);
        fd = shm_open(name, O_RDWR, S_IRUSR | S_IWUSR);
        ftruncate(fd, size);
        volatile char *p = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (p == MAP_FAILED)
        {
            printf("error map 1\n");
            return;
        }
        SHM_flags_Send[i] = (volatile struct SHM_flag_between_source_dest *)(p + intra_rank * sizeof(struct SHM_flag_between_source_dest));
        // if(intra_rank == 15)
        //     printf("%d -> %d \n",((int*)(SHM_flags_Send[i]))[0],((int*)(SHM_flags_Send[i]))[1] );
    }
    MPI_Barrier(Comm_intra);
}

void mread(int source, char *recvbuf, int size)
{
    //i am 0
    for (int start = 0; start < size; start += bufsize)
    {
        int localsize = (bufsize < (size - start)) ? bufsize : (size - start);
        //等待接收条件
        volatile int *front = &(SHM_flags_Recv[source]->front);
        int tail = SHM_flags_Recv[source]->tail;
        while (tail == (*front))
            ; //等待条件接收

        //开始接收
        volatile char *p = SHM_bufs_Recv[source] + tail * bufsize;
        for (int j = 0; j < localsize; j++)
            recvbuf[start + j] = p[j];
        //printf("check recv %c\n",((char *)recvbuf)[start]);
        //完成接收
        //puts("check recv flag");
        SHM_flags_Recv[source]->tail = (tail + 1) % bufN_between_ab;
    }
    // printf("check recv");
}
void mwrite(int target, char *sendbuf, int size)
{
    //i am 1
    for (int start = 0; start < size; start += bufsize)
    {
        int localsize = (bufsize < (size - start)) ? bufsize : (size - start);
        //等待发送条件
        int front = SHM_flags_Send[target]->front;            //front 只有写入的人可以改变，故而可以读到寄存器中
        volatile int *tail = &(SHM_flags_Send[target]->tail); //消费者会意外改变tail的值，故而用volatile指针
        while ((front + 1) % bufN_between_ab == (*tail))
            ;
        //进行发送
        volatile char *p = SHM_bufs_Send[target] + front * bufsize;
        for (int j = 0; j < localsize; j++)
            p[j] = sendbuf[start + j];
        //printf("check send %c\n",p[0]);
        front = (front + 1) % bufN_between_ab;
        //SHM_flags_Send[target]->size[front] = localsize;
        SHM_flags_Send[target]->front = front;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, procn;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procn);
    intra_procn = 32;
    intra_rank = rank % intra_procn;
    Comm_intra = MPI_COMM_WORLD;

    numa_noden = 1 ;
    //printf("numa_noden =%d\n",numa_noden);
    int a, b;
    a = atoi(argv[1]);
    b = atoi(argv[2]);
    unsigned int count = 0;
    GLEX_set_affinity(intra_rank);

    char *buf;
    int buf_size = 1;
    int ret = posix_memalign((void **)&buf, 4096, buf_size);
    if (ret)
    {
        fprintf(stderr, "posix_memalign: %s\n",
                strerror(ret));
        return -1;
    }
    const unsigned long nodeMask = 1UL << numa_node_of_cpu(intra_rank);
    mbind((void *)buf, buf_size, MPOL_BIND, &nodeMask, sizeof(nodeMask), MPOL_MF_MOVE | MPOL_MF_STRICT);
    memset(buf, 0, buf_size);
    int loopN = 10000000;
    // Init_SHMBUF();
    // Init_SHEMFLAG();

    // {
    //     double startT = MPI_Wtime();
    //     for (int loop = 0; loop < loopN; loop++)
    //     {
    //         if (rank == a)
    //         {
    //             //char buf[] = "ABCDEFGHIJKLMN";
    //             mwrite(b, buf, buf_size);
    //         }
    //         else if (rank == b)
    //         {
    //             //char buf[2];
    //             mread(a, buf, buf_size);
    //         }
    //     }
    //     double endT = MPI_Wtime();
    //     double totalT = (endT - startT) / loopN;
    //     double re = 0.0;
    //     MPI_Reduce(&totalT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, Comm_intra);
    //     if (intra_rank == 0)
    //         printf("%f us\t", re);
    // }

    for (int source = 0; source < 32; source++)
    {
        for (int target = 0; target < 32; target++)
        {
            if (source == target)
            {
                if (intra_rank == 0)
                    printf("0.0\t");
                continue;
            }
            if (rank == a)
            {
                //char buf[] = "ABCDEFGHIJKLMN";
                MPI_Send(buf, buf_size, MPI_CHAR, b, 0, Comm_intra);
            }
            else if (rank == b)
            {
                //char buf[2];
                MPI_Status status;
                MPI_Recv(buf, buf_size, MPI_CHAR, a, 0, Comm_intra, &status);
            }
            //开始测试MPI通信开销
            {
                double startT = MPI_Wtime();
                MPI_Barrier(Comm_intra);
                for (int loop = 0; loop < loopN; loop++)
                {
                    if (rank == a)
                    {
                        //char buf[] = "ABCDEFGHIJKLMN";
                        MPI_Send(buf, buf_size, MPI_CHAR, b, 0, Comm_intra);
                    }
                    else if (rank == b)
                    {
                        //char buf[2];
                        MPI_Status status;
                        MPI_Recv(buf, buf_size, MPI_CHAR, a, 0, Comm_intra, &status);
                    }
                }
                MPI_Barrier(Comm_intra);
                double endT = MPI_Wtime();
                double totalT = (endT - startT) * 1e6 / loopN;
                double re = 0.0;
                MPI_Reduce(&totalT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, Comm_intra);
                if (intra_rank == 0)
                    printf("%f \t", re);
            }
        }
        if (intra_rank == 0)
            puts("");
    }

    //puts("check");
    MPI_Finalize();
    return 0;
}