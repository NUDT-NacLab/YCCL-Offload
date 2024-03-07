
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

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    GLEXCOLL_Init(argc, argv);
    int rank, procn;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procn);
    int a, b, c, d;
    int size = (1 << 27);
    char *buf = (char *)malloc(size);
    for (size_t i = 0; i < size; i++)
    {
        buf[i] = 0;
    }
    int loopN = 500;
    {
        for (int a = 0; a < 32; a++)
        {
            for (int b = 0; b < 32; b++)
            {
                if (a == b)
                    continue;
                MPI_Status status;
                double time_start1 = MPI_Wtime();
                MPI_Barrier(MPI_COMM_WORLD);
                for (int loop = 0; loop < loopN; loop++)
                {
                    MPI_Barrier(MPI_COMM_WORLD);
                    if (rank == a)
                        MPI_Send(buf, size, MPI_CHAR, b, 0, MPI_COMM_WORLD);
                    if (rank == b)
                        MPI_Recv(buf, size, MPI_CHAR, a, 0, MPI_COMM_WORLD, &status);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                double time_end1 = MPI_Wtime();
                double time_local1 = (time_end1 - time_start1) / (loopN);
                double time_global1 = 0.0;
                MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                if (rank == 0)
                {
                    printf("%d %d %f\n", a, b, size / (1e6 * time_global1));
                    fflush(stdout);
                }
            }
        }
    }
    //测试节点内两条消息的竞争
    if (0)
    {
        for (int a = 0; a < 32; a += 4)
        {
            for (int b = 0; b < 32; b += 4)
            {
                if (a == b)
                    continue;
                for (int c = 1; c < 32; c += 4)
                {
                    for (int d = 1; d < 32; d += 4)
                    {
                        if (c == d)
                            continue;

                        double time_start1 = MPI_Wtime();
                        MPI_Barrier(MPI_COMM_WORLD);
                        for (int i = 0; i < loopN; i++)
                        {
                            MPI_Request req1, req2, req3, req4;
                            MPI_Status status;
                            if (rank == a)
                            {
                                MPI_Isend(buf, size, MPI_CHAR, b, 0, MPI_COMM_WORLD, &req1);
                            }
                            if (rank == b)
                            {
                                MPI_Irecv(buf, size, MPI_CHAR, a, 0, MPI_COMM_WORLD, &req2);
                            }
                            if (rank == c)
                            {
                                MPI_Isend(buf, size, MPI_CHAR, d, 0, MPI_COMM_WORLD, &req3);
                            }
                            if (rank == d)
                            {
                                MPI_Irecv(buf, size, MPI_CHAR, c, 0, MPI_COMM_WORLD, &req4);
                            }
                            if (rank == a)
                                MPI_Wait(&req1, &status);
                            if (rank == b)
                                MPI_Wait(&req2, &status);
                            if (rank == c)
                                MPI_Wait(&req3, &status);
                            if (rank == d)
                                MPI_Wait(&req4, &status);
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        double time_end1 = MPI_Wtime();
                        double time_local1 = (time_end1 - time_start1) / (loopN);
                        double time_global1 = 0.0;
                        MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                        if (rank == 0)
                        {
                            printf("%d %d %d %d %f\n", a, b, c, d, size / (1e6 * time_global1));
                            fflush(stdout);
                        }
                    }
                }
            }
        }
    }
    //     a = 4*atoi(argv[1]);
    //     b = 4*atoi(argv[2]);
    //     c = 4*atoi(argv[3])+1;
    //     d = 4*atoi(argv[4])+1;
    //     if(a == b || c==d)
    //     {
    //         MPI_Finalize();
    //         return 0;
    //     }
    //     unsigned loopN = 100000000;
    //     unsigned int count = 0;

    //     if(rank == 0)
    //     {
    //         GLEX_set_affinity(a);
    //         int fd = shm_open("GLEXshm2",O_CREAT | O_RDWR,0666);
    //         ftruncate(fd,64);
    //         char *p = (char *)mmap(NULL,64,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    //         memset(p,'S',64);
    //         MPI_Barrier(MPI_COMM_WORLD);
    //         //开始数据传输
    //         double t= 0;
    //         double start = MPI_Wtime();
    //         while(*p != 'R')
    //         {
    //             count++;
    //             t += ((SHM_FLIT *)p)->payload.T.data_double[0];
    //             *p = 'S';
    //             if(count == loopN)
    //                 break;
    //         }
    //         double end = MPI_Wtime();
    //         MPI_Status  status;
    //         double re1 = 0.0;
    //         MPI_Recv(&re1,1,MPI_DOUBLE,2,0,MPI_COMM_WORLD,&status);
    //         printf("%lf \t %lf\n",loopN/(end - start),re1);
    //         munmap(p,64);
    //     }
    //     else if(rank == 1)
    //     {
    //         GLEX_set_affinity(b);
    //         int fd = shm_open("GLEXshm2",O_CREAT | O_RDWR,0666);
    //         ftruncate(fd,64);
    //         char *p = (char *)mmap(NULL,64,PROT_READ,MAP_SHARED,fd,0);
    //         MPI_Barrier(MPI_COMM_WORLD);
    //         //开始数据传输
    //         while(*p != 'S')
    //         {
    //             ((SHM_FLIT *)p)->payload.size=1;
    //             ((SHM_FLIT *)p)->payload.T.data_double[0]=0.1;
    //             *p = 'R';
    //             count++;
    //             if(count == loopN)
    //                 break;
    //         }
    //         munmap(p,64);
    //     }
    //     else if(rank == 2)
    //     {
    //         GLEX_set_affinity(c);
    //         int fd = shm_open("GLEXshm1",O_CREAT | O_RDWR,0666);
    //         ftruncate(fd,64);
    //         char *p = (char *)mmap(NULL,64,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    //         memset(p,'S',64);
    //         MPI_Barrier(MPI_COMM_WORLD);
    //         //开始数据传输
    //         double t= 0;
    //         double start = MPI_Wtime();
    //         while(*p != 'R')
    //         {
    //             t += ((SHM_FLIT *)p)->payload.T.data_double[0];
    //             *p = 'S';
    //             count++;
    //             if(count == loopN)
    //                 break;
    //         }
    //         double end = MPI_Wtime();
    //        //printf("%lf \t msgs/second\n",loopN/(end - start));
    //         double re1 = loopN/(end - start);
    //         MPI_Send(&re1,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    //         munmap(p,64);
    //     }
    //     else if(rank == 3)
    //     {
    //         GLEX_set_affinity(d);
    //         int fd = shm_open("GLEXshm1",O_CREAT | O_RDWR,0666);
    //         ftruncate(fd,64);
    //         char *p = (char *)mmap(NULL,64,PROT_READ,MAP_SHARED,fd,0);
    //         MPI_Barrier(MPI_COMM_WORLD);
    //         //开始数据传输
    //         while(*p != 'S')
    //         {
    //             ((SHM_FLIT *)p)->payload.size=1;
    //             ((SHM_FLIT *)p)->payload.T.data_double[0]=0.1;
    //             *p = 'R';
    //             count++;
    //             if(count == loopN)
    //                 break;

    //         }
    //         munmap(p,64);
    //     }
    //   else{
    //         MPI_Barrier(MPI_COMM_WORLD);
    //     }
    MPI_Barrier(MPI_COMM_WORLD);
    GLEXCOLL_Finalize();
    MPI_Finalize();
    return 0;
}
