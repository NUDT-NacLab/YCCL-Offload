
#include <random>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>

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

double fun(int calc)
{
    if (calc == 0)
    {
        return 0.0;
    }
    else
    {
        int n = 1 << calc;
        int i = 0;
        double var = 0.0;
        for (i = 0; i < n; i += 1)
        {
            if (i % 2 == 0)
                var = sin(var);
            else
                var = cos(var);
        }
        return var;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int locR;
    MPI_Comm_rank(MPI_COMM_WORLD,&locR);
    MPI_Barrier(MPI_COMM_WORLD);
    if(locR == 0)
        fprintf(stderr,"================================lets start===================================\n");
    MPI_Barrier(MPI_COMM_WORLD);

    int calc = atoi(argv[1]);
    int sizeV = 1 << atoi(argv[2]);
    if (argc >= 4)
        CorePerNuma = 1 << atoi(argv[3]);
    else
    {
        CorePerNuma = 16;
    }
    if (CorePerNuma <= 0 || CorePerNuma > 32)
    {
        puts("CorePerNuma value error");
        exit(0);
    }
    if (argc >= 5)
        Childn_K = atoi(argv[4]);
    else
        Childn_K = 2;

    //if (Childn_K <= 0 || Childn_K > 15)
    //{
    //    puts("Childn_K value error");
    //    exit(0);
    //}
    GLEXCOLL_Init(argc, argv);
    MPI_Barrier(MPI_COMM_WORLD);
    //if (global_rank == 0)
    //    puts("check 78");
    if (ppn > 0)
        GLEXCOLL_InitAllreduce();
    MPI_Barrier(MPI_COMM_WORLD);
    //if (global_rank == 0)
    //    puts("check 74");
    // MPI_Finalize();
    // return EXIT_SUCCESS;
    // printf("CorePerNuma =%d\n ",CorePerNuma);
    // puts("GLEXCOLL_Init finished");
    // Get the size of the communicator
    int size = 0;
    MPI_Request req;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // if(size != 4)
    // {
    //     printf("This application is meant to be run with 4 MPI processes.\n");
    //     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // }
    // Get my rank
    int my_rank;
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    srand(my_rank);
    // if(my_rank == 0)
    //     printf("CorePerNuma = %d\n",CorePerNuma);
    double *source = (double *)malloc(sizeof(double) * sizeV);
    double *source1 = (double *)malloc(sizeof(double) * sizeV);

    // printf("sizeof(Payload) = %d\n",sizeof(struct Payload)); 64

    // int sources[sizeV];
    //随机数产生，
    static default_random_engine e;
    uniform_real_distribution<double> u0(-0.01, 0.01);
    double range = u0(e);
    uniform_real_distribution<double> u1(-range, range);
    uniform_int_distribution<unsigned> uinit(0, 1000);

    int i;
    //int loopN_GLEX = 10;
    int loopN_GLEX = 30000;
    int loopN_MPI = 30000;
    double *targets = (double *)malloc(sizeof(*targets) * sizeV);
    double *targets1 = (double *)malloc(sizeof(*targets1) * sizeV);

    // int targets[sizeV];
    //  Each MPI process sends its rank to reduction, root MPI process collects the result
    // warmup
    cout.precision(30); //设置精度
    int warmupN = 10;

    MPI_Barrier(MPI_COMM_WORLD);
    // GLEXCOLL_Iallreduce(source1, targets1, sizeV,MPI_DOUBLE,MPI_SUM);
    // GLEXCOLL_Wait_Iallreduce();
    // puts("check 122");
    // GLEXCOLL_AllreduceFinalize();
    // GLEXCOLL_Finalize();
    // MPI_Finalize();

    MPI_Barrier(MPI_COMM_WORLD);
    // if (global_rank == 0)
    // printf("rank=%d \n", global_rank);
    // if ((testmode & 2) == 2)
    {
        //测试MPI_IAllreduce
        double time_local1 = 0.0;
        for (int i = 0; i < loopN_MPI + 100; i++)
        {
            double time_start1 = MPI_Wtime();
            MPI_Allreduce(source, targets, sizeV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (i >= 100)
                time_local1 += (MPI_Wtime() - time_start1);
            MPI_Barrier(MPI_COMM_WORLD);
            // MPI_Iallreduce(source, targets, sizeV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req);
            // fun(calc);
            // MPI_Wait(&req, &status);
        }
        time_local1 = time_local1 / loopN_MPI;
        double time_global1 = 0.0;
        MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        time_global1 /= global_procn;
        if (my_rank == 0)
        {
            printf("MPI: %f\t\n", time_global1 * 1e6);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // if ((testmode & 1) == 1)
    {
        // for(INTRA_ALLREDUCE_TYPE = 3;INTRA_ALLREDUCE_TYPE<=6;INTRA_ALLREDUCE_TYPE+=3)
        {

            inter_comm_mode = OFFLOAD; //测试non offload,NON_OFFLOAD
            MPI_Barrier(MPI_COMM_WORLD);
            // GLEXCOLL_Iallreduce
	     srand(time(NULL));
	     for (int j = 0; j < sizeV; j++)
             {
	          source[j] = source1[j] = rand()/(RAND_MAX+1.0);
	     }
            for (int treeid = 0; treeid < 1 /*TreeNumber*/; treeid++)
            {
                _TreeID = 0;
                inter_comm_mode = 1;
                double time_local1 = 0.0;
                for (int i = 0; i < loopN_GLEX + 100; i++)
                //for (int i = 0; i < loopN_GLEX ; i++)
                {
	/*		for (int j = 0; j < sizeV; j++)
                     {
                        source[j] = source1[j] = rand()/(RAND_MAX+1.0);
                     }
                        usleep(uinit(e));*/
                    double time_start1 = MPI_Wtime();
                    GLEXCOLL_Iallreduce(source1, targets1, sizeV, MPI_DOUBLE, MPI_SUM);
                    GLEXCOLL_Wait_Iallreduce();
                    if (global_rank == -1)
		    {
		    	printf("%.30f\n",*targets1);
		    }
		    if (i >= 100)
                        time_local1 += (MPI_Wtime() - time_start1);
                    MPI_Barrier(MPI_COMM_WORLD);
                    // double recheck = global_procn * (i % 2);
                    // for (int j = 0; j < sizeV; j++)
                    // {
                    //     if (targets1[j] - recheck > 0.0001 || recheck - targets1[j] > 0.00001)
                    //     {
                    //         printf("结果验证错误：recheck = %f re=%f\n", recheck, targets1[j]);
                    //         exit(0);
                    //     }
                    // }
                    // MPI_Barrier(MPI_COMM_WORLD);
                    // if(i % 10000 == 0)
                }
                double time_end = MPI_Wtime();

                double time_local = time_local1 / loopN_GLEX;
                double time_global = 0.0;
                MPI_Reduce(&time_local, &time_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                time_global /= global_procn;
                if (my_rank == 0)
                {
                    printf("glexColl: %f\t\n", time_global * 1e6);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    // puts("");
    // puts("check 245");
    MPI_Barrier(MPI_COMM_WORLD);
    // puts("check 245");

    GLEXCOLL_AllreduceFinalize();
    GLEXCOLL_Finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
