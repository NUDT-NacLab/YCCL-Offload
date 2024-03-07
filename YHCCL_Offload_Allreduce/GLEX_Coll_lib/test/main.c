#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>
#include "glexcoll.h"

double fun(int calc)
{
    if(calc == 0)
    {
        return 0.0;
    }else{
        int n = 1<<calc;
        int i =0;
        double var = 0.0;
        for(i =0;i<n;i+=1)
        {
            var = sin(var);
            var = cos(var);
        }
        return var;
    }
}
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    //puts("lets start");
    int calc    = atoi(argv[1]);
    int sizeV = 1 << atoi(argv[2]);
    if(argc >= 4)
        CorePerNuma = 1<<atoi(argv[3]);
    else
    {
        CorePerNuma = 16;
    }
    if(CorePerNuma <= 0 || CorePerNuma>32)
    {
        puts("CorePerNuma value error");
        exit(0);
    }
    if(argc >= 5)
        Childn_K = atoi(argv[4]);
    else
        Childn_K = 2;
    
    if(Childn_K <= 0 || Childn_K > 15)
    {
        puts("CorePerNuma value error");
        exit(0);
    }
    GLEXCOLL_Init(argc, argv);
    //printf("CorePerNuma =%d\n ",CorePerNuma);
    //puts("GLEXCOLL_Init finished");
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
    // if(my_rank == 0)
    //     printf("CorePerNuma = %d\n",CorePerNuma);
    int *sources = (int *)malloc(sizeof(int) * sizeV);
    int *sources1 = (int *)malloc(sizeof(int) * sizeV);
    //int sources[sizeV];
    for(int i = 0;i<sizeV;i++)
        sources[i] =sources1[i] = 1;//4 + my_rank*my_rank;
    int *targets = (int *)malloc(sizeof(*targets) * sizeV);
    int *targets1 = (int *)malloc(sizeof(*targets1) * sizeV);
    //int targets[sizeV];
    int warmupN = 100;
    int loopN = 1000;
    
    // Each MPI process sends its rank to reduction, root MPI process collects the result
    //warmup
    for(int loop =0;loop<warmupN;loop++)
    {
        // if(my_rank == 0)
        //     printf("loop  =%d\n",loop);
        GLEXCOLL_Iallreduce(sources1, targets1, sizeV,MPI_INT,MPI_SUM);
        //puts("GLEXCOLL_Iallreduce posted");
        GLEXCOLL_Wait_Iallreduce();
        MPI_Iallreduce(sources,targets,sizeV,MPI_INT,MPI_SUM,MPI_COMM_WORLD,&req);
        MPI_Wait(&req,&status);
        if(1)
       {
            //正确性检查
            int re = 0;
            for(int i = 0;i<nprocs;i++)
            {
                re+=1;//(4 + i*i);
            }
            if(re != targets1[0])
            {
                printf("GLEXCOLL:loop = %d my_rank = %d\t ××××××××××× result %d(reduce) != %d(calc) \n",loop,my_rank,targets1[0],re);
                fflush(stdout);
            }else{
                // printf("my_rank = %d\t √ √ √ √ √ √ result %d(reduce) == %d(calc) \n",my_rank,targets1[0],re);
                // fflush(stdout);
            }
            if(re != targets[0])
            {
                printf("MPI:my_rank = %d\t ××××××××××× result %d(reduce) != %d(calc) \n",my_rank,targets1[0],re);
                fflush(stdout);
            }else{

            }
       }
    }
    
    // if(my_rank == 0)
    //     puts("check result");
    
    
    //测试MPI_IAllreduce
    double time_start1 = MPI_Wtime();
    for (int i = 0; i < loopN; i++)
    {
        MPI_Iallreduce(sources, targets, sizeV, MPI_INT, MPI_SUM, MPI_COMM_WORLD,&req);
        // Do some other job
        fun(calc);
        MPI_Wait(&req,&status);
    }
    double time_end1 = MPI_Wtime();

    double time_local1 = (time_end1 - time_start1) / loopN;
    double time_global1 = 0.0;
    MPI_Reduce(&time_local1, &time_global1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    //GLEXCOLL_Iallreduce
    double time_start = MPI_Wtime();
    for (int i = 0; i < loopN; i++)
    {
        GLEXCOLL_Iallreduce(sources1, targets1, sizeV,MPI_INT,MPI_SUM);
        fun(calc);
        GLEXCOLL_Wait_Iallreduce();
    }
    double time_end = MPI_Wtime();

    double time_local = (time_end - time_start) / loopN;
    double time_global = 0.0;
    MPI_Reduce(&time_local, &time_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    


    if(my_rank == 0)
    {
        printf("%f\t",time_global1*1e6);
        printf("%f\n",time_global*1e6);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    GLEXCOLL_Finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
