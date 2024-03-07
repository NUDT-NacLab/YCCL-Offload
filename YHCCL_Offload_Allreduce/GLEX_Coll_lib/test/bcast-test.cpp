
#include <random>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>

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

int main(int argc, char *argv[])
{
    int size_start = 0;
    int size_end = 26;
    int size_max = 1 + (1 << size_end);
    MPI_Init(&argc, &argv);
    Childn_K = (atoi(argv[2]));
    GLEXCOLL_Init(argc, argv);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    CorePerNuma = 1;
    _TreeID = 0;
    //Bcast 和 Allreduce共用初始化接口
    //测试每节点单进程
    GLEXCOLL_InitAllreduce();
    extern int SoftWare_Allreduce_Algorithm;
    extern int allreduce_slice_num;
    //puts("lets start");
    int argv1 = (atoi(argv[1]));
    double *sendbuf_MPI = (double *)malloc(size_max * sizeof(double));
    double *sendbuf_GLEX = (double *)malloc(size_max * sizeof(double));

    if (global_rank == 0)
    {
        for (int i = 0; i < size_max; i++)
        {
            sendbuf_GLEX[i] = (global_rank + 1.0) * i + 200;
        }
    }

    for (int size = size_start; size <= size_end; size++)
    {
        int count = 1 + (1 << size);
        int loopN = 0;
        if (size < 10)
            loopN = 2000;
        else if(size < 18)
            loopN = 400;
        else 
            loopN = 20;
            if (my_rank == 0)
                printf("double数组长度=%d \t", count);
       

        {
            //测试MPI
            double reT = 0.0;
            MPI_Barrier(MPI_COMM_WORLD);
            double startT = MPI_Wtime();
            for (int loop = 0; loop < loopN; loop++)
            {
                // MPI_Barrier(MPI_COMM_WORLD);
                MPI_Bcast(sendbuf_MPI, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            double endT = MPI_Wtime();
            // MPI_Barrier(MPI_COMM_WORLD);
            reT += (endT - startT);
            reT /= loopN;
            double re = 0.0;
            MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (my_rank == 0)
                printf("MPI:%.6f ", count, re);
        }
        for(int round=1;round<2;round++)
        {
            extern int Bcast_algorithm;
            if(round == 0)
	            Bcast_algorithm = K_ary_broadcast;//K_ary_broadcast_pipeline;
            else if(round == 1)
	            Bcast_algorithm = K_ary_broadcast_pipeline;
            //测试GLEX_Bcast
            double reT = 0.0;
            MPI_Barrier(MPI_COMM_WORLD);
                double startT = MPI_Wtime();
            for (int loop = 0; loop < loopN; loop++)
            {

                GLEXCOLL_Bcast(sendbuf_GLEX, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            double endT = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);
            reT += (endT - startT);
            reT /= loopN;
            double re = 0.0;
            MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (my_rank == 0)
                if(round == 0)
                    printf("K_ary_broadcast:%.6f ", re);
                else if(round == 1)
                    printf("K_ary_broadcast_pipeline:%.6f ", re);
            // puts("check 100");
            {
                //正确性验证
                {
#pragma parallel for
                    for (int i = 0; i < count; i++)
                    {
                        // printf("%f %f\n",sendbuf_MPI[i],sendbuf_GLEX[i]);
                        if (abs(i + 200.0 - sendbuf_GLEX[i]) > 0.001)
                        {
                            printf("广播结果错误，rank=%d index=%d\n", global_rank,i);
                            exit(0);
                        }
                    }
                }
            }
            // puts("check 116");
        }



        

        if (my_rank == 0)
            puts("");
    }

    free(sendbuf_MPI);
    free(sendbuf_GLEX);
    GLEXCOLL_AllreduceFinalize();
    MPI_Barrier(MPI_COMM_WORLD);
    GLEXCOLL_Finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}