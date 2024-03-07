
#include <mpi.h>
extern "C"
{
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
}

char *ptr;
int size = 64;
int loopN = 2000000;
int numa_noden;
int read1()
{
    int fd;
    //分配缓冲区，其地址为64字节对齐
    volatile double *bufRecv;
    posix_memalign((void **)&bufRecv, 64, 64 * loopN);
    const unsigned long nodeMask =-1;// 1UL << numa_node_of_cpu(_GLEXCOLL.intra_rank);
    mbind((void *)bufRecv, 64 * loopN, MPOL_BIND, &nodeMask, sizeof(nodeMask), MPOL_MF_MOVE | MPOL_MF_STRICT);
    for (int i = 0; i < 8 * loopN; i++)
        bufRecv[i] = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    fd = shm_open("GLEXSHM", O_RDWR, S_IRUSR | S_IWUSR);
    if (fd < 0)
    {
        printf("error open region\n");
        return 0;
    }
    ftruncate(fd, size);
    ptr = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (ptr == MAP_FAILED)
    {
        printf("error map\n");
        return 0;
    }
    int count = 0;
    double start = MPI_Wtime();
    while (count < loopN)
    {
        while (*ptr != 'R')
            ;
        for (int i = 0; i < 7; i++)
            bufRecv[(count << 3) + i] = ((SHM_FLIT *)ptr)->payload.T.data_double[i];
        *ptr = 'S';
        ++count;
    }
    free((void *)bufRecv);
    double end = MPI_Wtime();
    double time1 = loopN / (end - start);
    printf("%lf msg/second\n", time1);
    return 0;
}
int write1()
{
    //分配缓冲区，其地址为64字节对齐
    volatile double *bufSend;
    posix_memalign((void **)&bufSend, 64, 64 * loopN);
    const unsigned long nodeMask = 1UL << numa_node_of_cpu(_GLEXCOLL.intra_rank);
    mbind((void *)bufSend, 64 * loopN, MPOL_BIND, &nodeMask, sizeof(nodeMask), MPOL_MF_MOVE | MPOL_MF_STRICT);
    for (int i = 0; i < 8 * loopN; i++)
        bufSend[i] = 0.0;

    int fd;
    fd = shm_open("GLEXSHM", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
    if (fd < 0)
    {
        printf("error open region\n");
        return 0;
    }
    ftruncate(fd, size);
    ptr = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (ptr == MAP_FAILED)
    {
        printf("error map\n");
        return 0;
    }
    memset(ptr, 'S', size);
    MPI_Barrier(MPI_COMM_WORLD);

    int count = 0;
    while (count < loopN)
    {
        while (*ptr != 'S')
            ;
        for (int i = 0; i < 7; i++)
            ((SHM_FLIT *)ptr)->payload.T.data_double[i] = bufSend[(count << 3) + i];
        ptr[0] = 'R';
        ++count;
    }
    free((void *)bufSend);
    return 0;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    GLEXCOLL_Init(argc, argv);
    int rank, procn;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procn);
    numa_noden = 1 + numa_max_node();
    //printf("numa_noden =%d\n",numa_noden);
    int a, b;
    a = atoi(argv[1]);
    b = atoi(argv[2]);
    unsigned int count = 0;
    // if(rank == 0)
    // {
    //     GLEX_set_affinity(a);
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     //write1();
    // }
    // else if(rank == 1)
    // {
    //     GLEX_set_affinity(b);
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     //read1();
    // }else{
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    //puts("check");
    //第一步：测试MPI 点对点结果,由于MPI中采用了缓冲区合并消息传输，故而，本实现的延迟高于MPI延迟
    // {
    //     double time_start1 = MPI_Wtime();
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     for(int loop = 0;loop <= loopN;++loop)
    //     {
    //         double x;
    //         MPI_Status status;
    //         if(rank == 0)
    //         {
    //             MPI_Send(&x,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
    //             MPI_Recv(&x,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&status);
    //         }else if(rank == 1)
    //         {
    //             MPI_Recv(&x,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
    //             MPI_Send(&x,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    //         }
    //         MPI_Barrier(MPI_COMM_WORLD);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     double time_end1 = MPI_Wtime();
    //     double time_local1 = (time_end1 - time_start1) / (2*loopN);
    //     double time_global1 = 0.0;
    //     MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //    if(rank == 0)
    //     {
    //         printf("%f\t",time_global1*1e6);
    //         fflush(stdout);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    //第二步：测试SDOT 点对点结果
    {
        double time_start1 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        for (int loop = 0; loop <= loopN; ++loop)
        {
            double x;
            if (rank == a)
            {
                volatile char *p = SHM_bufs_Send[b];
                while (*p != 'S')
                    ;
                ((volatile SHM_FLIT *)p)->payload.T.data_double[0] = x;
                *p = 'R';
                while (*p != 'B')
                    ;
                x = ((volatile SHM_FLIT *)p)->payload.T.data_double[0];
                *p = 'S';
            }
            else if (rank == b)
            {
                volatile char *p = SHM_bufs_Recv[a];
                while (*p != 'R')
                    ;
                x = ((volatile SHM_FLIT *)p)->payload.T.data_double[0];
                ((volatile SHM_FLIT *)p)->payload.T.data_double[0] = x;
                *p = 'B';
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double time_end1 = MPI_Wtime();
        double time_local1 = (time_end1 - time_start1) / (2 * loopN);
        double time_global1 = 0.0;
        MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            printf("%f \t", time_global1 * 1e9);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (0)
    {
        double time_start1 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        for (int loop = 0; loop <= loopN; ++loop)
        {
            double x;
            if (rank == a)
            {
                volatile char *p = SHM_bufs_Send[b];
                volatile char *p1 = SHM_bufs_Recv[b];
                while (*p != 'S')
                    ;
                ((volatile SHM_FLIT *)p1)->payload.T.data_double[0] = x;
                *p = 'R';
                while (*p != 'B')
                    ;
                x = ((volatile SHM_FLIT *)p1)->payload.T.data_double[0];
                *p = 'S';
            }
            else if (rank == b)
            {
                volatile char *p = SHM_bufs_Recv[a];
                volatile char *p1 = SHM_bufs_Send[a];
                while (*p != 'R')
                    ;
                x = ((volatile SHM_FLIT *)p1)->payload.T.data_double[0];
                ((volatile SHM_FLIT *)p1)->payload.T.data_double[0] = x;
                *p = 'B';
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double time_end1 = MPI_Wtime();
        double time_local1 = (time_end1 - time_start1) / (2 * loopN);
        double time_global1 = 0.0;
        MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            printf("%f \t", time_global1 * 1e9);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (0)
    {
        double time_start1 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        for (int loop = 0; loop <= loopN; ++loop)
        {
            double x;
            //ping
            if (rank == a)
            {
                volatile char *p = SHM_bufs_Send[b];
                while (*p != 'S')
                    ;
                ((volatile SHM_FLIT *)p)->payload.T.data_double[0] = x;
                *p = 'R';
            }
            else if (rank == b)
            {
                volatile char *p = SHM_bufs_Recv[a];
                while (*p != 'R')
                    ;
                x = ((volatile SHM_FLIT *)p)->payload.T.data_double[0];
                *p = 'S';
            }
            //pong
            if (rank == a)
            {
                volatile char *p = SHM_bufs_Recv[b];
                while (*p != 'R')
                    ;
                x = ((volatile SHM_FLIT *)p)->payload.T.data_double[0];
                *p = 'S';
            }
            else if (rank == b)
            {
                volatile char *p = SHM_bufs_Send[a];
                while (*p != 'S')
                    ;
                ((volatile SHM_FLIT *)p)->payload.T.data_double[0] = x;
                *p = 'R';
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double time_end1 = MPI_Wtime();
        double time_local1 = (time_end1 - time_start1) / (2 * loopN);
        double time_global1 = 0.0;
        MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            printf("%f\t", time_global1 * 1e9);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    {
        double time_start1 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        for (int loop = 0; loop <= loopN; ++loop)
        {
            double x;
            //ping
            if (rank == a)
            {
                volatile char *p = SHM_bufs_Send[b];
                volatile char *p1 = SHM_bufs_Recv[b];
                while (*p != 'S')
                    ;
                __sync_synchronize();
                ((volatile SHM_FLIT *)p1)->payload.T.data_double[0] = x;
                __sync_synchronize();
                *p = 'R';
            }
            else if (rank == b)
            {
                volatile char *p = SHM_bufs_Recv[a];
                volatile char *p1 = SHM_bufs_Send[a];
                while (*p != 'R')
                    ;
                __sync_synchronize();
                x = ((volatile SHM_FLIT *)p1)->payload.T.data_double[0];
                __sync_synchronize();
                *p = 'S';
            }
            // //pong
            if (rank == a)
            {
                volatile char *p = SHM_bufs_Recv[b];
                volatile char *p1 = SHM_bufs_Send[b];
                while (*p != 'R')
                    ;
                __sync_synchronize();
                x = ((volatile SHM_FLIT *)p1)->payload.T.data_double[0];
                __sync_synchronize();
                *p = 'S';
            }
            else if (rank == b)
            {
                volatile char *p = SHM_bufs_Send[a];
                volatile char *p1 = SHM_bufs_Recv[a];
                while (*p != 'S')
                    ;
                __sync_synchronize();
                ((volatile SHM_FLIT *)p1)->payload.T.data_double[0] = x;
                __sync_synchronize();
                *p = 'R';
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double time_end1 = MPI_Wtime();
        double time_local1 = (time_end1 - time_start1) / (2 * loopN);
        double time_global1 = 0.0;
        MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            printf("%f\t", time_global1 * 1e9);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    {
        double time_start1 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        static int pvt_release_state = 0;
        static int pvt_gather_state = 0;
        for (int loop = 0; loop <= loopN; ++loop)
        {
            double x;
            //ping
            pvt_release_state++;
            pvt_gather_state++;
            if (rank == a)
            {
                *(shared_release_flags[a]) = pvt_release_state;
                __sync_synchronize();
                while (*(shared_gather_flags[b]) != pvt_gather_state)
                    ;
                __sync_synchronize();
                x = ((volatile double *)Reduce_buffers[b])[0];
            }
            else if (rank == b)
            {
                while (*(shared_release_flags[a]) != pvt_release_state)
                    ;
                __sync_synchronize();
                ((volatile double *)Reduce_buffers[b])[0] = x;
                __sync_synchronize();
                *(shared_gather_flags[b]) = pvt_gather_state;
            }
            //pong
            pvt_release_state++;
            pvt_gather_state++;
            if (rank == a)
            {
                while (*(shared_release_flags[b]) != pvt_release_state)
                    ;
                __sync_synchronize();
                ((volatile double *)Reduce_buffers[a])[0] = x;
                __sync_synchronize();
                *(shared_gather_flags[a]) = pvt_gather_state;
            }
            else if (rank == b)
            {
                *(shared_release_flags[b]) = pvt_release_state;
                __sync_synchronize();
                while (*(shared_gather_flags[a]) != pvt_gather_state)
                    ;
                __sync_synchronize();
                x = ((volatile double *)Reduce_buffers[a])[0];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double time_end1 = MPI_Wtime();
        double time_local1 = (time_end1 - time_start1) / (2 * loopN);
        double time_global1 = 0.0;
        MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            printf("%f\n", time_global1 * 1e9);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // if(rank == 0)
    // {
    //
    //     int fd = shm_open("GLEXshm",O_CREAT | O_RDWR,S_IRUSR | S_IWUSR);
    //     ftruncate(fd,64);
    //     char *p = (char *)mmap(NULL,64,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    //     memset(p,'S',64);
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     //开始数据传输
    //     double t= 0;
    //     double start = MPI_Wtime();
    //     while(count != loopN){
    //         while(*p != 'R') ;

    //         {
    //             count++;
    //             t += ((SHM_FLIT *)p)->payload.T.data_double[0];
    //             *p = 'S';
    //         }
    //     }
    //     double end = MPI_Wtime();
    //     printf("%lf \t msgs/second\n",loopN/(end - start));
    //     munmap(p,64);
    // }
    // else if(rank == 1)
    // {
    //
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     int fd = shm_open("GLEXshm",O_CREAT | O_RDWR,S_IRUSR | S_IWUSR);
    //     ftruncate(fd,64);
    //     char *p = (char *)mmap(NULL,64,PROT_READ,MAP_SHARED,fd,0);
    //     if   (p == MAP_FAILED) {
    //         printf ( "error map\n" );
    //         return   0;
    //     }
    //     //开始数据传输
    //     while(count != loopN)
    //     {
    //         while(*p != 'S');
    //         puts("loop 00");
    //         {
    //             p[63] = 'B';
    //             puts("loop0");
    //             ((SHM_FLIT *)p)->payload.T.data_double[0]=0.1;
    //             *p = 'R';
    //             count++;
    //         }
    //     }
    //     munmap(p,64);

    // }else{
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    GLEXCOLL_Finalize();
    MPI_Finalize();
    return 0;
}
