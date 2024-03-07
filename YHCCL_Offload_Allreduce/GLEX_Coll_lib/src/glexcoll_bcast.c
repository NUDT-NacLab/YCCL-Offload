#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <endian.h>
#include <sched.h>
#include <errno.h>
#include <mpi.h>
#include <stdint.h>
#include <omp.h>
#include <sys/mman.h>
#include <pthread.h>
#include "glex.h"
#include "glexcoll.h"
#include "glexallreduce.h"

//逻辑rank只支持每节点单进程
static int my_logical_rank;
static int group_rank;
static int group_procn;
static int group_intra_rank;
static int group_intra_procn;
static int group_inter_rank;
static int group_inter_procn;
int Bcast_algorithm;
int slice_size = (1<<16);
//根据逻辑rank的值获取comm中真实rank的值
inline int logical_to_real_rank(int logical_rank, int procn, int root)
{
	return (logical_rank + procn + root) % procn;
}
//根据真实rank值获取逻辑rank值
inline int real_to_logical_rank(int real_rank, int procn, int root)
{
	return (real_rank + procn - root) % procn;
}
//获取K-ary数中i号孩子。
inline int Get_logical_childi(int logical_rank, int procn, int i, int childK)
{
	return (childK * logical_rank + 1 + i);
}
//获取K-ary完全树中孩子的起始逻辑rank
inline int Get_logical_child_start(int logical_rank, int childK)
{
	return (logical_rank * childK + 1);
}

//获取K-ary完全树中孩子的结束逻辑rank
inline int Get_logical_child_end(int logical_rank, int childK)
{
	return ((logical_rank + 1) * childK);
}
//获取K-ary完全树中孩子的父节点逻辑rank
inline int Get_logical_parent(int logical_rank, int childK)
{
	if (logical_rank == 0)
		return 0;
	else
		return (logical_rank - 1) / childK;
}

void glexcoll_K_ary_broadcast(void *data_p, int size, int real_root, MPI_Comm comm)
{
	static int flag_count = 0;
	MPI_Request reqv[32];
	MPI_Status statusV[32];
	if (my_logical_rank != 0)
	{
		// printf("Childn_K=%d\n",Childn_K);
		//除了root以外的节点都需要接收消息
		int logical_parent = Get_logical_parent(my_logical_rank, Childn_K);
		int real_parent = logical_to_real_rank(logical_parent, group_procn, real_root);
		MPI_Recv(data_p, size, MPI_CHAR, real_parent, flag_count, comm, statusV);
		// int recvn;
		// MPI_Get_count(statusV,MPI_CHAR,&recvn);
		// if(recvn != size)
		// {
		// 	puts("recv check error");
		// 	exit(0);
		// }
		// printf("%d's parent=%d\n", global_rank, real_parent);
	}
	int child_logical_start = Get_logical_child_start(my_logical_rank, Childn_K);
	if (child_logical_start < group_procn)
	{
		//代表当前进程还有孩子存在
		int child_logical_end = mmin(group_procn - 1, Get_logical_child_end(my_logical_rank, Childn_K));
		int count = 0;
		for (int logical_child = child_logical_start; logical_child <= child_logical_end; logical_child++)
		{
			int real_child = logical_to_real_rank(logical_child, group_procn, real_root);
			// printf("%d's logical child=%d real child=%d\n", global_rank, logical_child, real_child);
			MPI_Isend(data_p, size, MPI_CHAR, real_child, flag_count, comm, &(reqv[count]));
			count++;
		}
		MPI_Waitall(count, reqv, statusV);
	}
	flag_count++;
	// MPI_Barrier(comm);
	// exit(0);
}
void glexcoll_K_ary_broadcast_pipeline(void *data_p, int size, int real_root, MPI_Comm comm)
{
	// puts("check");
	extern int mmin(int a, int b);
	//用于发送和接受使用的缓冲区
	void *buf_send = data_p;
	void *buf_recv = data_p;
	static MPI_Request reqv[32];
	static MPI_Status statusV[32];
	//计算该消息的传输需要多少次迭代
	int roundn = size / slice_size;
	if (size % slice_size != 0)
		roundn++;
	//第一步需要获取自身的节点类型：ROOT，MID，LEAF。
	int my_type = 0;
	int child_logical_start = Get_logical_child_start(my_logical_rank, Childn_K);
	if (my_logical_rank == 0)
	{
		my_type = 0; //ROOT
	}
	else if (child_logical_start < group_procn)
	{
		my_type = 1; //MID节点
	}
	else
	{
		my_type = 2; //LEAF节点
	}
	//第二步：对于非root节点，进行消息接收
	int remain_size = size;
	int round_size = mmin(remain_size, slice_size);
	remain_size -= round_size;

	int logical_parent = Get_logical_parent(my_logical_rank, Childn_K);
	int real_parent = logical_to_real_rank(logical_parent, group_procn, real_root);

	if (my_type == 1 || my_type == 2)
	{
		//MID和LEAF需要接收
		MPI_Recv(buf_recv, round_size, MPI_CHAR, real_parent, 0, comm, statusV);
		buf_recv += round_size;
	}
	while (remain_size > 0)
	{
		int count_req = 0;
		if (my_type != 2)
		{
			//对于root和MID需要将消息发送出去。
			int child_logical_end = mmin(group_procn - 1, Get_logical_child_end(my_logical_rank, Childn_K));
			for (int logical_child = child_logical_start; logical_child <= child_logical_end; logical_child++)
			{
				int real_child = logical_to_real_rank(logical_child, group_procn, real_root);
				// printf("%d's logical child=%d real child=%d\n", global_rank, logical_child, real_child);
				MPI_Isend(buf_send, round_size, MPI_CHAR, real_child, 0, comm, &(reqv[count_req++]));
				// MPI_Isend(data_p, (1<<27), MPI_DOUBLE, real_child, 0, comm, &(reqv[count_req]));
			}
			buf_send += round_size;
		}
		{
			//重新计算下一阶段的round_size;
			round_size = mmin(remain_size, slice_size);
			remain_size -= round_size;
		}
		if (my_type != 0)
		{
			//对于MID和LEAF来说。需要接受来自parent的消息。
			MPI_Irecv(buf_recv, round_size, MPI_CHAR, real_parent, 0, comm, &(reqv[count_req++]));
			buf_recv += round_size;
		}
		MPI_Waitall(count_req, reqv, statusV);
	}
	if (my_type != 2)
	{
		int count_req = 0;
		//对于root和MID需要将消息发送出去。
		int child_logical_end = mmin(group_procn - 1, Get_logical_child_end(my_logical_rank, Childn_K));
		for (int logical_child = child_logical_start; logical_child <= child_logical_end; logical_child++)
		{
			int real_child = logical_to_real_rank(logical_child, group_procn, real_root);
			// printf("%d's logical child=%d real child=%d\n", global_rank, logical_child, real_child);
			MPI_Isend(buf_send, round_size, MPI_CHAR, real_child, 0, comm, &(reqv[count_req]));
			count_req++;
			// MPI_Isend(data_p, (1<<27), MPI_DOUBLE, real_child, 0, comm, &(reqv[count_req]));
		}
		buf_send += round_size;
		MPI_Waitall(count_req, reqv, statusV);
	}
	// exit(0);
}

//广播函数的要求必须是comm中的进程数整除节点数
void GLEXCOLL_Bcast(
	void *data_p,
	int count,
	MPI_Datatype datatype,
	int source_proc,
	MPI_Comm comm)
{
	MPI_Comm_rank(comm, &group_rank);
	MPI_Comm_size(comm, &group_procn);
	group_intra_procn = ppn;
	group_intra_rank = group_rank % ppn;
	group_inter_rank = group_rank / ppn;
	group_inter_procn = group_procn / ppn;

	my_logical_rank = real_to_logical_rank(group_rank, group_procn, source_proc);
	int msg_size = 8;
	// MPI_Type_size(datatype, &msg_size);
	msg_size *= count;
	//算法选择器
	switch (Bcast_algorithm)
	{
		case K_ary_broadcast:
			/* code */
			glexcoll_K_ary_broadcast(data_p, msg_size, source_proc, comm);
			break;
		case K_ary_broadcast_pipeline:
			glexcoll_K_ary_broadcast_pipeline(data_p, msg_size, source_proc, comm);
			break;
		default:
			MPI_Bcast(data_p, count, datatype, source_proc, comm);
			break;
	}
}