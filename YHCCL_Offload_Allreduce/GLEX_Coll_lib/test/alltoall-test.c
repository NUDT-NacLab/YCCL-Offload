#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "glexcoll.h"
#include "glexalltoall.h"

int main(int argc, char *argv[])
{
    int procn, rank;
    int iter_end = 9;
    int iter_start = 1;
    int block_size = 1 << iter_end; //block_size<=10
    int loopN = 1000000;

    double *send_buff;
    double *recv_buff;
    double *send_buff1;
    double *recv_buff1;
    extern int leaderN;
    extern int num_of_ongoing_msg;
    extern int RT_Self_Adapting_on;
    extern int Alltoall_algorithm; //BRUCK;
    extern int Memcpy_On;
    extern int Num_of_ongoing_and_SelfAdapting_ON;
    leaderN = atoi(argv[1]);
    Memcpy_On = 1;
    int num_of_ongoing_msg_max = (atoi(argv[2]));
    num_of_ongoing_msg = num_of_ongoing_msg_max;
    MPI_Init(&argc, &argv);
    GLEXCOLL_Init(argc, argv);

    glexcoll_InitAlltoall();
    // puts("check 21");

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procn);

    send_buff = (double *)malloc(block_size * procn * sizeof(double));
    recv_buff = (double *)malloc(block_size * procn * sizeof(double));
    send_buff1 = send_buff; //(double *)malloc(block_size * procn * sizeof(double));
    recv_buff1 = (double *)malloc(block_size * procn * sizeof(double));

    for (int size = iter_start; size <= iter_end; size++)
    {
        //size=17;
        block_size = 1 << size;
        long long total_size = procn * procn * block_size;
        if (size < 7)
            loopN = 1000;
        else if (size < 13)
            loopN = 10;
        for (int i = 0; i < procn; i++)
        {
            for (int j = 0; j < block_size; j++)
            {
                send_buff1[i * block_size + j] = send_buff[i * block_size + j] = rank * 100000.0 + i;
            }
        }
        if (global_rank == 0)
            printf("%d * double\t", block_size);

        //warmUP mpi
        {
            MPI_Barrier(MPI_COMM_WORLD);
            double startT = MPI_Wtime();
            extern int mmin(int a, int b);
            for (int loop = 0; loop < mmin(loopN, 10); loop++)
            {
                MPI_Alltoall(send_buff, block_size, MPI_DOUBLE, recv_buff, block_size, MPI_DOUBLE, MPI_COMM_WORLD);
                //
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        //测试MPI
        {
            double reT = 0.0;
            MPI_Barrier(MPI_COMM_WORLD);
            for (int loop = 0; loop < loopN; loop++)
            {
                double startT = MPI_Wtime();
                MPI_Alltoall(send_buff, block_size, MPI_DOUBLE, recv_buff, block_size, MPI_DOUBLE, MPI_COMM_WORLD);
                double endT = MPI_Wtime();
                reT += (endT - startT);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            reT /= loopN;
            double re = 0.0;
            MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (rank == 0)
                printf("%f ", re, block_size);
            {
                for (int proc = 0; proc < procn; proc++)
                {
                    for (int m = 0; m < block_size; m++)
                    {
                        double val = proc * 100000.0 + rank;
                        if (recv_buff[proc * block_size + m] - val > 0.01 || recv_buff[proc * block_size + m] - val < -0.01)
                        {
                            printf("测试MPI:check error rank%d \n", rank);
                            exit(0);
                        }
                    }
                }
            }
        }
        //测试glexcoll_alltoall
        extern int COMMUNICATION_ON;
        extern int Intra_Gather_Scatter_ON;
        extern int MATRIX_Tanfform_ON;

        extern int slice_size;
        slice_size = (1 << 16);
        COMMUNICATION_ON = 1;
        Intra_Gather_Scatter_ON = 1;
        MATRIX_Tanfform_ON = 1;
        //Warmup
        //puts("start warmup");
        // if (loopN != 1)
        // {
        //     num_of_ongoing_msg = 1;
        //     RT_Self_Adapting_on = 1;
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     for (int loop = 0; loop < 10; loop++)
        //     {
        //         GLEXCOLL_Alltoall(send_buff1, block_size * sizeof(double), recv_buff1, block_size * sizeof(double), MPI_COMM_WORLD);
        //         MPI_Barrier(MPI_COMM_WORLD);
        //     }
        // }
        // puts("finish warmup");
        //for(num_of_ongoing_msg = 1;num_of_ongoing_msg <=6;num_of_ongoing_msg++)
        Memcpy_On = 1;
        //Alltoall_algorithm = DIRECT;//DIRECT_Kleader_NODE_AWARE_RDMA;
        //for(int num_of_ongoing_msg = 1;num_of_ongoing_msg<=4;num_of_ongoing_msg*=2)
        num_of_ongoing_msg = 2;
        {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            //slice_size =(1<<c);
            double reT = 0.0;
            for (int loop = 0; loop < loopN; loop++)
            {
                double startT = MPI_Wtime();
                GLEXCOLL_Alltoall(send_buff1, block_size * sizeof(double), recv_buff1, block_size * sizeof(double), MPI_COMM_WORLD, MPI_DOUBLE);
                double endT = MPI_Wtime();
                reT += (endT - startT);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            reT = reT / loopN;
            double re = 0.0;
            MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            // if (rank == 0)
            //     printf("%lf %ld GB %lf GB/s\n", re,8*total_size/(1e9),8*total_size/(re*1e9));

            if (rank == 0)
                printf("%lf\t", re);
            {
                for (int proc = 0; proc < procn; proc++)
                {
                    for (int m = 0; m < block_size; m++)
                    {
                        double val = proc * 100000.0 + rank;
                        if (recv_buff1[proc * block_size + m] - val > 0.01 || recv_buff1[proc * block_size + m] - val < -0.01)
                        {
                            printf("测试glexcoll_alltoall:check error rank%d %f!=%f \n", rank,recv_buff1[proc * block_size + m],val);
                            exit(0);
                        }
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        // for(int leadn=1;leadn<ppn;leadn++)
        // {

        // }
        //测试节点感知算法
        if (0)
            for (int t = 0; t < 1; ++t)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                double startT = MPI_Wtime();
                for (int loop = 0; loop < loopN; loop++)
                {
                    GLEXCOLL_Alltoall(send_buff1, block_size * sizeof(double),
                                      recv_buff1, block_size * sizeof(double), MPI_COMM_WORLD, MPI_DOUBLE);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                double endT = MPI_Wtime();
                double reT = (endT - startT) / loopN;
                double re = 0.0;
                MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                if (rank == 0)
                    printf("%f \t", re);
                {
                    for (int proc = 0; proc < procn; proc++)
                    {
                        for (int m = 0; m < block_size; m++)
                        {
                            double val = proc * 100000.0 + rank;
                            if (recv_buff1[proc * block_size + m] - val > 0.01 || recv_buff1[proc * block_size + m] - val < -0.01)
                            {
                                // printf("测试glexcoll_alltoall:check error rank%d \n", rank);
                                // exit(0);
                            }
                        }
                    }
                }
            }

        if (rank == 0)
            puts("");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(send_buff);
    free(recv_buff);
    GLEXCOLL_AlltoallFinalize();
    GLEXCOLL_Finalize();
    MPI_Finalize();
    return 0;
}