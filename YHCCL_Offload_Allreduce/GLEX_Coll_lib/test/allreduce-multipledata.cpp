
#include <random>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <unistd.h>
extern "C"
{
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "glexcoll.h"
#include "glexallreduce.h"
}
using namespace std;
extern int leaderN;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    int iphrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &iphrank);
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, 30000);
    generator.seed(iphrank);

    leaderN = 1;
    if (iphrank == 0)
        puts("check 24");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    GLEXCOLL_Init(argc, argv);
    MPI_Barrier(MPI_COMM_WORLD);
    if (global_rank == 0)
        puts("check 27");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    CorePerNuma = 1;
    Childn_K = (atoi(argv[2]));
    _TreeID = 0;
    GLEXCOLL_InitAllreduce();
    MPI_Barrier(MPI_COMM_WORLD);
    if (global_rank == 0)
        puts("check 32");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    extern int allreduce_slice_num;
    int argv1 = (atoi(argv[1]));
    int size_max = 1 << argv1;
    double *sendbuf_MPI = (double *)malloc(size_max * sizeof(double));
    double *recvbuf_MPI = (double *)malloc(size_max * sizeof(double));
    double *sendbuf_GLEX = (double *)malloc(size_max * sizeof(double));
    double *recvbuf_GLEX = (double *)malloc(size_max * sizeof(double));

    static default_random_engine e;
    uniform_real_distribution<double> u0(-10, 10);
    double range = u0(e);
    for (int i = 0; i < argv1; i++)
    {
        sendbuf_MPI[i] = sendbuf_GLEX[i] = 1.0; // u0(e);
    }
    int loopN;
    int loopN_MPI;
    // puts("check");
    for (int sz = 0; sz <= argv1; sz++)
    {

        if (sz <= 8)
            loopN = 100000;
        else
            loopN = 6000;
        loopN_MPI = loopN;
        // loopN = 2;
        MPI_Barrier(MPI_COMM_WORLD);
        int size = 1 << sz;
        {
            //测试mpi
            double reT = 0.0;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 0; i < loopN_MPI; i++)
            {
                double startT = MPI_Wtime();
                MPI_Allreduce(sendbuf_MPI, recvbuf_MPI, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                double endT = MPI_Wtime();
                reT += (endT - startT);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                // if(global_rank == 0) puts("check 113");
                // MPI_Reduce(sendbuf_MPI, recvbuf_MPI, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            reT /= loopN_MPI;
            double re = 0.0;
            MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (global_rank == 0)
                printf("消息大小：%d \t MPI_Allreduce: %f us\t", size, re * 1e6);
            fflush(stdout);
        }
        // {
        //     //测试mpi
        //     double reT = 0.0;
        //     for (int i = 0; i < loopN; i++)
        //     {
        //         double startT = MPI_Wtime();
        //         // MPI_Allreduce(sendbuf_MPI, recvbuf_MPI, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //         MPI_Reduce(sendbuf_MPI, recvbuf_MPI, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //         double endT = MPI_Wtime();
        //         reT += (endT - startT);
        //         MPI_Barrier(MPI_COMM_WORLD);
        //         MPI_Barrier(MPI_COMM_WORLD);
        //     }
        //     reT /= loopN;
        //     double re = 0.0;
        //     MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        //     if (global_rank == 0)
        //         printf("1<<%d\t %f \t", sz, re * 1e6);
        // }

        SoftWare_Allreduce_Algorithm = Small_message_node_aware_POSIX_RDMA;
        // Small_message_node_aware_POSIX_RDMA;
        //  Recursize_doubling_OMP; //Recursize_doubling_OMP  ; // Recursize_doubling_OMP; //K_nomial_tree_OMP //Recursize_doubling_Slicing
        extern int vec_self_tuned_two_dimentional_allreduce_comm_NUM;
        extern int self_tuned_tow_dimentional_allreduce_DimX;
        // loopN = 1;
        if (0)
            for (self_tuned_tow_dimentional_allreduce_DimX = 1; self_tuned_tow_dimentional_allreduce_DimX <= inter_procn / 2; self_tuned_tow_dimentional_allreduce_DimX++)
            {
                {
                    //测试GLEX
                    double reT = 0.0;
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    for (int i = 0; i < loopN; i++)
                    {
                        double startT = MPI_Wtime();
                        GLEXCOLL_Allreduce(sendbuf_GLEX, recvbuf_GLEX, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        double endT = MPI_Wtime();
                        reT += (endT - startT);
                        MPI_Barrier(MPI_COMM_WORLD);
                        // if(global_rank == 0) puts("check 113");
                        // GLEXCOLL_Allreduce_NonOffload(sendbuf_GLEX, size, recvbuf_GLEX);
                    }
                    reT /= loopN;
                    double re = 0.0;
                    MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                    if (global_rank == 0)
                        printf("DimX = %d %f \n", self_tuned_tow_dimentional_allreduce_DimX, re * 1e6);
                    fflush(stdout);
                }
                {
//正确性检查
#pragma omp parallel for
                    for (int i = 0; i < size; i++)
                    {
                        if (abs(recvbuf_GLEX[i] - recvbuf_MPI[i]) > 0.001)
                        {
                            printf("glex result: %f != MPI result: %f\n", recvbuf_GLEX[i], recvbuf_MPI[i]);
                            fflush(stdout);

                            exit(0);
                        }
                    }
                }
            }
        MPI_Barrier(MPI_COMM_WORLD);
        // SoftWare_Allreduce_Algorithm = K_nomial_tree_OMP;
        // for(allreduce_slice_num=(1<<10);allreduce_slice_num<=(1<<size_max);allreduce_slice_num*=2)
        extern int allreduce_k;
        // for(int t= 10;t<=argv1;t++)
        //  for (allreduce_k = 2; allreduce_k <= 64; allreduce_k *= 2)
        allreduce_k = 4;
        {
            // allreduce_slice_num = (1<<t);
            MPI_Barrier(MPI_COMM_WORLD);
            {
                //测试GLEX
                double reT = 0.0;
                MPI_Barrier(MPI_COMM_WORLD);
                // warmup
                // if (0)
                for (int x = 1; x < 400; x++)
                {
                    for (int i = 0; i < size; i++)
                        ((volatile double *)sendbuf_GLEX)[i] = x & i;
                    // if (intra_rank == 1)
                    //     usleep(distribution(generator));
                    GLEXCOLL_Allreduce(sendbuf_GLEX, recvbuf_GLEX, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    {
                        //正确性检查
                        for (int i = 0; i < size; i++)
                        {
                            if (abs(recvbuf_GLEX[i] - (x & i) * global_procn) > 0.001)
                            {
                                // if(intra_rank ==0)
                                printf("size=%d 正确性检查未通过 x=%d global rank  %d i=%d re = %f\n", size, x, global_rank, i, recvbuf_GLEX[i]);
                                fflush(stdout);
                                exit(0);
                            }
                        }
                    }

                    // if (global_rank == distribution(generator) % global_procn)
                    //     usleep(distribution(generator));
                    // MPI_Barrier(MPI_COMM_WORLD);

                    // if (global_rank == 0)
                    // {
                    //     printf(" check loop = %d\n", x);
                    // }
                    // MPI_Barrier(MPI_COMM_WORLD);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                double startT = MPI_Wtime();
                for (int i = 0; i < loopN; i++)
                {
                    // *(volatile double*)recvbuf_GLEX = 0.0;
                    GLEXCOLL_Allreduce(sendbuf_GLEX, recvbuf_GLEX, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    // *(volatile double*)sendbuf_GLEX = 0.0;

                    // MPI_Barrier(MPI_COMM_WORLD);
                    // MPI_Barrier(MPI_COMM_WORLD);
                }
                double endT = MPI_Wtime();
                reT += (endT - startT);
                reT /= loopN;
                double re = 0.0;
                MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                if (global_rank == 0)
                    printf(" GLEX_Allreduce: %f us", re * 1e6);
                fflush(stdout);
            }
            if (0)
            {
                //正确性检查
                for (int i = 0; i < size; i++)
                {
                    if (abs(recvbuf_GLEX[i] - recvbuf_MPI[i]) > 0.001)
                    {
                        // if(intra_rank ==0)
                        printf("check error b global rank  %d re = %f\n", global_rank, recvbuf_GLEX[i]);
                        exit(0);
                    }
                }
            }
        }
        if (global_rank == 0)
            puts("");
        fflush(stdout);
    }

    free(sendbuf_MPI);
    free(recvbuf_MPI);
    free(sendbuf_GLEX);
    free(recvbuf_GLEX);
    GLEXCOLL_AllreduceFinalize();
    GLEXCOLL_Finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
