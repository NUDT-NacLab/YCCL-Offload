#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include "glexcoll.h"
#include "glexalltoall.h"
void intra_memory_barrier_alltoall();
#define DATA_TYPE int
int main(int argc, char *argv[])
{
	// puts("check 101");
	int procn, rank;
	int iter_end = 19;
	int iter_start = 0;

	int block_size = 1 << iter_end; //block_size<=10
	int loopN = 1000000;
	// puts("check");
	DATA_TYPE *send_buff;
	DATA_TYPE *recv_buff;
	DATA_TYPE *send_buff1;
	DATA_TYPE *recv_buff1;
	extern int leaderN;
	extern int num_of_ongoing_msg;
	extern int RT_Self_Adapting_on;
	extern int Alltoall_algorithm; //BRUCK;
	extern int Memcpy_On;
	extern int Num_of_ongoing_and_SelfAdapting_ON;
	leaderN = atoi(argv[1]);
	Memcpy_On = 1;
	int num_of_ongoing_msg_max = 2; //(atoi(argv[2]));
	num_of_ongoing_msg = num_of_ongoing_msg_max;
	//  puts("check 34");
	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	//  puts("check 31");
	// return 0;
	GLEXCOLL_Init(argc, argv);
	MPI_Barrier(MPI_COMM_WORLD);
	// puts("check 36");
	glexcoll_init_alltoall_shared_memory_buffer();
	int rowN = global_procn; //128;//(int)sqrt(global_procn);
	int color = (global_rank) / rowN + 1;
	if (global_rank == 0)
	{
		// printf("rowN = %d\n", rowN);
	}
	MPI_Comm row_commm;
	MPI_Comm_split(MPI_COMM_WORLD, color, global_rank, &row_commm);
	// puts("check 34");
	glexcoll_InitAlltoall_new(row_commm);
	MPI_Barrier(MPI_COMM_WORLD);
	// if (global_rank == 0)
	// 	puts("check 35");
	MPI_Comm_rank(row_commm, &rank);
	MPI_Comm_size(row_commm, &procn);
	for (int i = 0; i < 1000; i++)
	{
		intra_memory_barrier_alltoall();
	}
	send_buff1 = (DATA_TYPE *)get_is_senddata_buffer(intra_rank); //send_buff; //(DATA_TYPE *)malloc(block_size * procn * sizeof(DATA_TYPE));
	recv_buff1 = (DATA_TYPE *)get_is_recvdata_buffer(intra_rank); //(DATA_TYPE *)malloc(block_size * procn * sizeof(DATA_TYPE));

	struct GLEXCOLL_a2a_bufmh a2a_mh;
	// glexcoll_register_alltoall_buffer(send_buff1, recv_buff1, &a2a_mh);
	glexcoll_register_alltoall_buffer_new(send_buff1, recv_buff1, block_size * procn * sizeof(DATA_TYPE), MPI_COMM_WORLD, &a2a_mh);

	// if (global_rank == 0)
	// 	puts("check 64");
	send_buff = send_buff1; //(DATA_TYPE *)malloc(block_size * procn * sizeof(DATA_TYPE));
	recv_buff = recv_buff1; //(DATA_TYPE *)malloc(block_size * procn * sizeof(DATA_TYPE));
	if (iter_end < iter_start)
	{
		puts(" iter_end < inter_start error happend");
		exit(0);
	}
	for (int size = iter_start; size <= iter_end; size++)
	{
		//size=17;
		block_size = 1 << size;
		// printf("check %d %d\n",iter_start,iter_end);
		long long int total_size = procn * procn * block_size;
		// printf("check1 %d %d\n",iter_start,iter_end);
		if (size < 10)
			loopN = 8000;
		else
			loopN = 400;
		// loopN = 2;
		if (global_rank == 0)
			printf("num of single msg: %8.0d  ", block_size);
		// puts("exit 69");
		// exit(0);
		//测试MPI
		MPI_Barrier(MPI_COMM_WORLD);
		// if(global_rank == 0)
		//     puts("check 79");

		for (int i = 0; i < procn; i++)
		{
			for (int j = 0; j < block_size; j++)
			{
				// if(global_rank == 79 && (i * block_size + j)%10000 == 0)
				//     printf("%d\n",i * block_size + j);
				send_buff[i * block_size + j] = color * 3000000 + rank * 1000 + i;
			}
		}

		double min_re = 199999999.0;
		MPI_Barrier(MPI_COMM_WORLD);
		// if(global_rank == 0){
		//     puts("check 110");
		// 	fflush(stdout);
		// }
		// if(0)
		{
			double reT = 0.0;
			MPI_Barrier(MPI_COMM_WORLD);
			for (int loop = 0; loop < loopN; loop++)
			{

				// for (int i = 0; i < procn; i++)
				// {
				//     for (int j = 0; j < block_size; j++)
				//     {
				//         send_buff[i * block_size + j] = rank * 100000 + i + loop;
				//     }
				// }
				double startT = MPI_Wtime();
				MPI_Alltoall(send_buff1, block_size, MPI_INT, recv_buff1, block_size, MPI_INT, row_commm);
				double endT = MPI_Wtime();
				reT += (endT - startT);
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				// {
				//     for (int proc = 0; proc < procn; proc++)
				//     {
				//         for (int m = 0; m < block_size; m++)
				//         {
				//             int val = proc * 100000 + rank + loop;
				//             if (recv_buff[proc * block_size + m] - val != 0)
				//             {
				//                 printf("测试MPI:check error rank%d \n", rank);
				//                 exit(0);
				//             }
				//         }
				//     }
				// }
			}
			MPI_Barrier(row_commm);
			reT /= loopN;
			double re = 0.0;
			MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			if (global_rank == 0)
				printf("%.9f ", re);
			// if (re < min_re)
			min_re = re;

			if (global_rank == 0)
			{
				long long sz = (block_size * 4) * (global_procn - 1);
				sz *= global_procn;
				double gbsize = sz / (1000000000.0);
				double bw = gbsize / min_re;
				printf("MPI网络吞吐率: %3.5f GB/s  ", bw);
			}

			// if (global_rank == 0)
			// {
			//     long long sz = (block_size * 4) * (global_procn - 1);
			//     sz *= global_procn;
			//     double gbsize = sz/(1000000000.0);
			//     double bw = gbsize / re;
			//     printf("%f GB %f GB/s %f",gbsize, bw, re);
			// }
		}

		MPI_Barrier(MPI_COMM_WORLD);
		// if(global_rank == 0)
		//     puts("check 166");

		MPI_Barrier(MPI_COMM_WORLD);
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
		// if(global_rank == 0)
		//     puts("finish warmup");
		//for(num_of_ongoing_msg = 1;num_of_ongoing_msg <=6;num_of_ongoing_msg++)
		Memcpy_On = 1;
		//Alltoall_algorithm = DIRECT;//DIRECT_Kleader_NODE_AWARE_RDMA;
		Intra_Gather_Scatter_ON = 1;
		MATRIX_Tanfform_ON = 1;
		extern int Topology_aware_up;
		Topology_aware_up = 0;
		Intra_Gather_Scatter_ON = 1;
		MATRIX_Tanfform_ON = 1;
		for (num_of_ongoing_msg = 1; num_of_ongoing_msg <= 8; num_of_ongoing_msg *= 2)
			for (int pjt = 0; pjt < 1; ++pjt)
			{
				// for(Topology_aware_up = 0;Topology_aware_up<1;Topology_aware_up+=1)
				Topology_aware_up = pjt;
				fflush(stdout);
				// if(global_rank == 0) printf("pjt = %d\n",pjt);
				MPI_Barrier(MPI_COMM_WORLD);
				if (pjt == 0)
					Alltoall_algorithm = L_a2a;
				{
					if (pjt == 1)
					{
						Alltoall_algorithm = NMPML;
						leaderN = 2;
					}
					if (pjt == 2)
					{
						Alltoall_algorithm = NMPML;
						leaderN = 4;
					}
					if (pjt == 3)
					{
						Alltoall_algorithm = ONMPML;
						leaderN = 2;
					}
					if (pjt == 4)
					{
						Alltoall_algorithm = ONMPML;
						leaderN = 4;
					}
					if (pjt == 5)
					{
						Alltoall_algorithm = SONMPML;
						leaderN = 2;
					}
					if (pjt == 6)
					{
						Alltoall_algorithm = SONMPML;
						leaderN = 4;
					}
					//warmup+正确性验证
					// if(0)
					for (int loop = 0; loop < 2; loop++)
					{
						for (int i = 0; i < procn; i++)
						{
							for (int j = 0; j < block_size; j++)
							{
								send_buff1[i * block_size + j] = rank * 100000 + i + loop * 10000000;
							}
						}
						// *recv_buff1 = -10;
						MPI_Barrier(row_commm);
						//GLEXCOLL_Alltoall_new(send_buff1, block_size * sizeof(DATA_TYPE), recv_buff1, block_size * sizeof(DATA_TYPE), row_commm, MPI_INT);
						GLEXCOLL_Alltoall_pjt(&a2a_mh, block_size * sizeof(DATA_TYPE));
						// *send_buff1 = -11;
						MPI_Barrier(row_commm);
						MPI_Barrier(MPI_COMM_WORLD);
						// if(global_rank == 0) puts("check 191");
						// if(0)
						{
							for (int proc = 0; proc < procn; proc++)
							{
								for (int m = 0; m < block_size; m++)
								{
									int val = proc * 100000 + rank + loop * 10000000;
									if (recv_buff1[proc * block_size + m] - val != 0)
									{
										printf("测试GLEXCOLL_a2a:check error rank %d ->%d %d!=%d\n", proc, rank, recv_buff1[proc * block_size + m], val);
										exit(0);
									}
								}
							}
						}
					}

					//     Alltoall_algorithm = TPDS17_Cache_oblivious_intra_node;
					// if (pjt == 1)
					//     Alltoall_algorithm = shared_memory_direct_intra_node;
					// if (pjt == 2)
					//     Alltoall_algorithm = TPDS17_Cache_oblivious_intra_node_NUMA;

					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Barrier(MPI_COMM_WORLD);
					//slice_size =(1<<c);

					//时间测试
					double reT = 0.0;
					MPI_Barrier(row_commm);
					for (int loop = 0; loop < loopN; loop++)
					{

						double startT = MPI_Wtime();
						// MPI_Alltoall(send_buff, block_size, MPI_INT, recv_buff, block_size, MPI_INT, row_commm);
						// GLEXCOLL_Alltoall_new(send_buff1, block_size * sizeof(DATA_TYPE), recv_buff1, block_size * sizeof(DATA_TYPE), row_commm, MPI_INT);
						GLEXCOLL_Alltoall_pjt(&a2a_mh, block_size * sizeof(DATA_TYPE));
						double endT = MPI_Wtime();
						reT += (endT - startT);
						MPI_Barrier(row_commm);
						MPI_Barrier(row_commm);
						// if(rank == 0)
						// {
						//     puts("-------------------------------------------PJT check------------------------------------------------");
						// }
						// if(global_rank == 0) puts("check 255");
					}
					MPI_Barrier(MPI_COMM_WORLD);
					reT = reT / loopN;
					double re = 0.0;
					MPI_Reduce(&reT, &re, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
					extern int mmin(int a, int b);
					// if( re < min_re)
					min_re = re;
					// if (rank == 0)
					//     printf("%lf %ld GB %lf GB/s\n", re,8*total_size/(1e9),8*total_size/(re*1e9));

					// if (global_rank == 0)
					// 	printf(" %6.6f ", re);
					for (int h = 0; h < block_size * procn; h++)
					{
						recv_buff1[h] = 0;
					}

					MPI_Barrier(MPI_COMM_WORLD);
				}
				if (0)
					if (Alltoall_algorithm == MPI_LINEAR_EXCHANGE_COMPRESS)
					{
						extern double alltoall_compress_ratio_sum;
						//如果使用通信压缩算法，统计出平均压缩率
						//alltoall_compress_ratio_sum->alltoall_compress_ratio_avg;
						alltoall_compress_ratio_sum /= global_procn;
						double ratio_avg;
						MPI_Reduce(&alltoall_compress_ratio_sum, &ratio_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
						ratio_avg /= global_procn;
						if (global_rank == 0)
						{
							printf("%f ", ratio_avg);
						}
					}
				// if(0)
				if (global_rank == 0)
				{
					long long sz = (block_size * 4) * (global_procn - 1);
					sz *= global_procn;
					double gbsize = sz / (1000000000.0);
					double bw = gbsize / min_re;
					printf("网络吞吐率: %4.5f GB/s ", bw);
				}
			}
		if (global_rank == 0)
			puts("");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// puts("314");
	// GLEXCOLL_AlltoallFinalize();
	// GLEXCOLL_Alltoall_Finalize_new();
	MPI_Barrier(MPI_COMM_WORLD);
	// puts("316");
	// GLEXCOLL_Finalize();
	// puts("318");

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
