#include <iostream>
#include <random>
#include <algorithm>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include "glexcoll.h"
#include "glexalltoall.h"

int *sendCounts;
int *recvCounts;
int *sdisp;
int *rdisp;

uint64_t *sendCounts64;
uint64_t *recvCounts64;
uint64_t *sdisp64;
uint64_t *rdisp64;

using namespace std;
extern int alltoallV_ongoingmsgN;
int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    GLEXCOLL_Init(argc, argv);
    // if(inter_rank ==0)
    // puts("start GLEXCOLL_InitAlltoallV");
    GLEXCOLL_InitAlltoallV(MPI_COMM_WORLD);
    alltoallV_ongoingmsgN = atoi(argv[1]);
    MPI_Barrier(MPI_COMM_WORLD);
    // if(inter_rank ==0)
    // puts("finish init");
    //分配MPI参数空间
    int Max_size = 1 << 22;
    sendCounts = new int[inter_procn];
    recvCounts = new int[inter_procn];
    sdisp = new int[inter_procn];
    rdisp = new int[inter_procn];

    //分配GLEX参数空间
    sendCounts64 = new uint64_t[inter_procn];
    recvCounts64 = new uint64_t[inter_procn];
    sdisp64 = new uint64_t[inter_procn];
    rdisp64 = new uint64_t[inter_procn];

    default_random_engine e;
    uniform_int_distribution<unsigned> u(0, Max_size);
    for (int i = 0; i < inter_procn; i++)
        sendCounts64[i] = sendCounts[i] = 1 << 22; //u(e);//10+inter_rank; //1 << 22; //
    MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < inter_procn; i++)
        recvCounts64[i] = recvCounts[i];
    sdisp64[0] = sdisp[0] = 0;
    for (int i = 1; i < inter_procn; ++i)
    {
        sdisp[i] = sendCounts[i - 1] + sdisp[i - 1];
        sdisp64[i] = sendCounts64[i - 1] + sdisp64[i - 1];
    }
    rdisp64[0] = rdisp[0] = 0;
    for (int i = 1; i < inter_procn; ++i)
    {
        rdisp[i] = recvCounts[i - 1] + rdisp[i - 1];
        rdisp64[i] = recvCounts64[i - 1] + rdisp64[i - 1];
    }
    for (int i = 1; i < inter_procn; ++i)
    {
        if ((sendCounts64[i] != sendCounts[i]))
            puts("error 初始化出错 (sendCounts64[i] != sendCounts[i])");
        if (recvCounts64[i] != recvCounts[i])
            puts("error 初始化出错 (recvCounts64[i] != recvCounts[i])");
        if (sdisp64[i] != sdisp[i])
            puts("error 初始化出错 (sdisp64[i] != sdisp[i])");
        if (rdisp64[i] != rdisp[i])
            puts("error 初始化出错 (rdisp64[i] != rdisp[i])");
    }
    //分配发送和接收缓冲区
    uint64_t Sendsize = 0, Recvsize = 0;
    for (int i = 0; i < inter_procn; i++)
    {
        Sendsize += sendCounts[i];
        Recvsize += recvCounts[i];
    }
    char *sendBuf = (char *)new char[Sendsize];
    char *sendBuf1 = (char *)new char[Sendsize];
    char *recvBuf = (char *)new char[Recvsize];
    char *recvBuf1 = (char *)new char[Recvsize];

    int sendrecvBufferID = 0;
    {
        GLECOLL_Register_SendRecv_Pair(sendBuf1, Sendsize, recvBuf1, Recvsize, MPI_COMM_WORLD);
        GLEXCOLL_Alltoallv_regist_recv_disps(rdisp64);
    }

//int PairId = GLECOLL_Register_SendRecv_Pair(sendBuf1, Sendsize, recvBuf1, Recvsize, MPI_COMM_WORLD) ;

//发送和接收缓冲区的初始化
#pragma omp parallel for
    for (int i = 0; i < Sendsize; i++)
    {
        sendBuf[i] = sendBuf1[i] = (inter_rank + i) % 26 + 'a';
    }

    // if(inter_rank ==0)
    //     puts("start MPI alltoallv rdma");
    //开始执行循环测试
    int loopN = 5;
    double time_max = 0.0;
    double time_min = 99999.0;
    {

        double time = 0.0;
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < loopN; i++)
        {
            double startT = MPI_Wtime();
            MPI_Alltoallv(sendBuf, sendCounts, sdisp, MPI_CHAR, recvBuf, recvCounts, rdisp, MPI_CHAR, Comm_inter);
            double endT = MPI_Wtime();
            time += (endT - startT);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            // if(inter_rank ==0)
            //     printf("MPI_Alltoallv %d\n",i);
        }
        time /= loopN;
        MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (inter_rank == 0)
            cout << time_max << "\t";
    }
    if (inter_rank == 0)
        time_min = std::min(time_min, time_max);

    MPI_Barrier(MPI_COMM_WORLD);
    // if(inter_rank ==0)
    //     puts("start MPI LINEAR_SHIFT");
    extern int AlltoallvAlgorithm;
    for (alltoallV_ongoingmsgN = 1; alltoallV_ongoingmsgN <= 10; alltoallV_ongoingmsgN += 1)
    // alltoallV_ongoingmsgN = 2;
    {
        if (inter_rank == 0)
            cout << endl;
        {
            AlltoallvAlgorithm = LINEAR_SHIFT;
            //开始关于Linear Shift 步重排的研究
            double time = 0.0;
            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 0; i < loopN; i++)
            {
                double startT = MPI_Wtime();
                GLEXCOLL_Alltoallv(sendBuf1, sendCounts64, sdisp64, recvBuf1, recvCounts64, rdisp64, MPI_CHAR);
                double endT = MPI_Wtime();
                time += (endT - startT);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                // if(global_rank == 0)
                //     puts("check finish round1");
                //MPI_Alltoallv(sendBuf1,sendCounts,sdisp,MPI_CHAR, recvBuf1,recvCounts,rdisp,MPI_CHAR,Comm_inter);
            }
            time /= loopN;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (inter_rank == 0)
                time_min = std::min(time_min, time_max);
            if (inter_rank == 0)
                cout << time_max << "\t";
        }
        {
            //开始基于RDMA实现的MPI_Alltoall
            MPI_Barrier(MPI_COMM_WORLD);
            AlltoallvAlgorithm = RDMA_LINEAR_SHIFT;
            double time = 0.0;
            static int c = 0;

            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 0; i < loopN; i++)
            {
                double startT = MPI_Wtime();
                //GLEXCOLL_Alltoallv(sendBuf1, sendCounts64, sdisp64, recvBuf1, recvCounts64, rdisp64);
                GLEXCOLL_Alltollv_RDMA(sendrecvBufferID, sendBuf1, recvBuf1, sendCounts64, sdisp64);
                double endT = MPI_Wtime();
                time += (endT - startT);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                // if(global_rank == 0)
                //     puts("check finish round1");
            }
            time /= loopN;
            MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (inter_rank == 0)
                time_min = std::min(time_min, time_max);
            if (inter_rank == 0)
                cout << time_max << "\t";
        }
        // {
        //     // GLEXCOLL_Alltoallv_Performance_detection(sendBuf1, sendCounts64, sdisp64, recvBuf1, recvCounts64, rdisp64);
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     double time = 0.0;
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     for (int i = 0; i < loopN; i++)
        //     {
        //         double startT = MPI_Wtime();
        //         GLEXCOLL_Alltoallv_RDMA_REORDER(sendrecvBufferID, sendBuf1, recvBuf1, sendCounts64, sdisp64);
        //         double endT = MPI_Wtime();
        //         time += (endT - startT);
        //         MPI_Barrier(MPI_COMM_WORLD);
        //         MPI_Barrier(MPI_COMM_WORLD);
        //     }
        //     time /= loopN;
        //     MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        //     if (inter_rank == 0)
        //         time_min = std::min(time_min, time_max);
        //     if (inter_rank == 0)
        //         cout << time_max << "\t";
        // }
        // {
        //     //基于数据压缩的RDMA通信
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     AlltoallvAlgorithm = RDMA_LINEAR_SHIFT;
        //     //开始关于RDMA的性能优化
        //     double time = 0.0;
        //     static int c = 0;

        //     MPI_Barrier(MPI_COMM_WORLD);
        //     for (int i = 0; i < loopN; i++)
        //     {
        //         double startT = MPI_Wtime();
        //         //GLEXCOLL_Alltoallv(sendBuf1, sendCounts64, sdisp64, recvBuf1, recvCounts64, rdisp64);
        //         GLEXCOLL_Alltollv_RDMA(sendrecvBufferID, sendBuf1, recvBuf1, sendCounts64, sdisp64);
        //         double endT = MPI_Wtime();
        //         time += (endT - startT);
        //         MPI_Barrier(MPI_COMM_WORLD);
        //         MPI_Barrier(MPI_COMM_WORLD);
        //         // if(global_rank == 0)
        //         //     puts("check finish round1");
        //     }
        //     time /= loopN;
        //     MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        //     if (inter_rank == 0)
        //         time_min = std::min(time_min, time_max);
        //     if (inter_rank == 0)
        //         cout << time_max << "\t";
        // }
    }
    {
        //统计总共发送和接受的消息字节数
        uint64_t total_bytes;

        MPI_Reduce(&Sendsize, &total_bytes, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        if (inter_rank == 0)
            cout << total_bytes / 1e9 << " GB\t" << total_bytes / (1e9 * time_min) << " GB/S\t";
        MPI_Barrier(MPI_COMM_WORLD);
    }
    {
        //结果检查
        int start = 0;
#pragma omp parallel for
        for (int i = 0; i < inter_procn; i++)
        {
            if (recvCounts[i] != recvCounts64[i])
            {

                puts("check error (recvCounts[i] != recvCounts64[i]) x");
                exit(0);
            }
            for (int j = 0; j < recvCounts[i]; j++)
            {
                int index = start + j;
                if (recvBuf1[index] != recvBuf[index])
                {
                    printf("error 结果检查错误 send_rank=%d -> recv_rank=%d\n ", i, inter_rank);
                    exit(0);
                }
                //exit(0);
            }
            start += recvCounts[i];
        }
#pragma omp parallel for
        for (int i = 0; i < Recvsize; i++)
            if (recvBuf1[i] != recvBuf[i])
            {
                printf("error 结果检查错误  recv_rank=%d\n ", inter_rank);
                exit(0);
            }
    }

    GLEXCOLL_Finalize();
    MPI_Barrier(Comm_inter);
    MPI_Barrier(Comm_inter);
    // puts("start MPI_Finalize");
    // puts("end MPI_Finalize");
}