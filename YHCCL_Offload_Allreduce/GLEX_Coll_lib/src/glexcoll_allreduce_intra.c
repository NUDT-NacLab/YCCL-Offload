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
#include <sys/mman.h>
#include <pthread.h>
#include "glex.h"
#include "glexcoll.h"
#include "glexallreduce.h"

volatile char *SHM_bufs_Recv[64];
volatile char *SHM_bufs_Send[64];
volatile char *SHM_bufs_bcast;
volatile char *SHM_bufs_bcast1;
void Intel_release(int bcastflag, void *sendbuf, void *recvbuf, int count);
void Intel_gather(int reduceflag, void *sendbuf, void *recvbuf, int count);

void one_double_allreduce_intra_l2_cache_aware_intel(void *sendbuf, void *recvbuf, int count)
{
	//第一步为gather
	static int pvt_gather_state = 0;
	static int pvt_release_state = 0;
	int block_leader = (intra_rank >> 2) << 2;
	double *re = (double *)recvbuf;
	for (int i = 0; i < count; i++)
		re[i] = ((double *)sendbuf)[i];
	if (block_leader == intra_rank)
	{
		pvt_gather_state++;
		//第一步规约到block_leader
		// sleep(1);
		for (int c = 1; c <= 3; ++c)
		{
			//对每个孩子
			volatile int *p = shared_gather_flags[intra_rank + c];
			while (*p != pvt_gather_state)
				;
			// printf("*p =%d\n",*p);
			for (int i = 0; i < count; i++)
				re[i] += ((double *)(Reduce_buffers[intra_rank + c]))[i];
		}
		*(shared_gather_flags[intra_rank]) = pvt_gather_state;
		// printf("%d %f\n",intra_rank,re[0]);
		/////////////////////////////////////////////////////////////////
		// pvt_gather_state++;
		// pvt_release_state++;
		// if(intra_rank == 0)
		// {
		// 	//root
		// 	for(int c = 1;c<(intra_procn>>2);++c)
		// 	{
		// 		volatile int * p = shared_gather_flags[(c<<2)];
		// 		while (*p != pvt_gather_state);
		// 		//puts("check 296");
		// 		for(int i = 0;i<count;i++)
		// 			re[i] += ((double *)(Reduce_buffers[(c<<2)]))[i];
		// 	}
		// 	for (int i = 0; i < count; i++)
		// 		((double *)(broad_cast_buffer))[i] = re[i];
		// 	*(shared_release_flags[intra_rank]) = pvt_release_state;
		// }else{
		// 	//block leader
		// 	//向上规约
		// 	for (int i = 0; i < count; i++)
		// 	{
		// 		((double *)Reduce_buffers[intra_rank])[i] = re[i];
		// 	}
		// 	*(shared_gather_flags[intra_rank]) = pvt_gather_state;
		// 	//等待广播结果
		// 	while(*(shared_release_flags[0]) != pvt_release_state);
		// 		for(int c = 0;c<count;c++)
		// 		{
		// 			re[c] = ((double *)broad_cast_buffer)[c];
		// 		}

		// }
		// pvt_gather_state--;
		// pvt_release_state--;
		////////////////////////////////////////////////////////////////////////////////
		//将结果广播回其它进程
		pvt_release_state++;
		for (int i = 0; i < count; i++)
			((double *)(broad_cast_buffer))[i] = re[i];
		*(shared_release_flags[intra_rank]) = pvt_release_state;
	}
	else
	{
		pvt_gather_state++;
		for (int i = 0; i < count; i++)
		{
			((double *)Reduce_buffers[intra_rank])[i] = ((double *)sendbuf)[i];
		}
		*(shared_gather_flags[intra_rank]) = pvt_gather_state;

		//等待广播结果
		pvt_release_state++;
		while (*(shared_release_flags[block_leader]) != pvt_release_state)
			;
		for (int c = 0; c < count; c++)
		{
			((double *)recvbuf)[c] = ((double *)broad_cast_buffer)[c];
		}
		// printf("%d %f\n",intra_rank,((double *)recvbuf)[0]);
	}
	// exit(0);
}

void one_double_allreduce_intra_REDUCE_BCAST_Intel(void *sendbuf, void *recvbuf, int count)
{
	Intel_release(0, sendbuf, recvbuf, count);
	Intel_gather(1, sendbuf, recvbuf, count);
	Intel_release(1, recvbuf, recvbuf, count);
	Intel_gather(0, recvbuf, recvbuf, count);
	//第一步为gather
	// static int pvt_gather_state = 111;
	// static int pvt_release_state = 111;
	// double * re = (double *) recvbuf;
	// for(int i = 0;i<count;i++) re[i] = ((double *)sendbuf)[i];
	// if(intra_rank == 0)
	// {
	// 	//root
	// 	pvt_gather_state++;
	//     for(int i = 1;i<intra_procn;++i)
	//     {
	//         while(*(shared_gather_flags[i]) != pvt_gather_state);
	//             for(int c = 0;c<count;c++)
	//                 re[c] += ((double *)Reduce_buffers[i])[c];
	//     }
	// 	pvt_release_state++;
	//     for(int c = 0;c<count;c++)
	//     {
	//         ((double *)broad_cast_buffer)[c] = re[c];
	//     }
	// 	*shared_release_flags[0] = pvt_release_state;

	// }else{
	// 	//childs
	// 	pvt_gather_state++;
	//     for(int i = 0;i<count;i++)
	//         {
	//             ((double *)Reduce_buffers[intra_rank])[i] = ((double *)sendbuf)[i];
	//         }
	// 	 *(shared_gather_flags[intra_rank]) = pvt_gather_state;
	// 	pvt_release_state++;
	// 	while (*(shared_release_flags[0]) != pvt_release_state) ;
	//     for(int c = 0;c<count;c++)
	//     {
	//         re[c] = ((double *)broad_cast_buffer)[c];
	//     }
	// }
	// printf("%f \n",re[0]);
	//
}
void one_double_allreduce_intra_REDUCE_to_0(void *sendbuf, void *recvbuf, int count)
{
	// __sync_synchronize();
	// double * re = recvbuf;
	// if(_GLEXCOLL.intra_rank == 0)
	// {

	// 	volatile char *p1 = (_GLEXCOLL.SHM_0+((ppn-1)<<6));
	// 	*p1 = 'P';
	// 	for(int i = 0;i<count;++i) re[i] = ((double *)sendbuf)[i];
	// 	//sleep(1);
	// 	for(int i = 0;i<ppn -1 ;++i)
	// 	{
	// 		volatile char *p = _GLEXCOLL.SHM_buf+(i<<6);
	//     	while (*p != 'R');
	// 		// printf("check %d\n",_GLEXCOLL.intra_rank);
	// 		// fflush(stdout);
	// 			for(int i = 0;i<count;i++)
	// 				re[i] += ((volatile SHM_FLIT *)p)->payload.T.data_double[i];
	// 	}
	// }else{
	// 	volatile char *p = _GLEXCOLL.SHM_buf;
	// 	*p ='S';
	//     while (*p != 'S');
	//         //((SHM_FLIT *)p)->payload.size=1;
	// 		for(int i = 0;i<count;i++)
	//         	((volatile SHM_FLIT *)p)->payload.T.data_double[i]=((double *)sendbuf)[i];
	// 		*p = 'R';
	// }
	// puts("check 187");
	int block_leader = (intra_rank >> 2) << 2;
	double *re = (double *)recvbuf;
	for (int i = 0; i < count; i++)
		re[i] = ((double *)sendbuf)[i];

	if (block_leader == intra_rank)
	{
		//第一步规约到block_leader
		for (int c = 1; c <= 3; ++c)
		{
			//对每个孩子
			volatile char *p = SHM_bufs_Recv[intra_rank + c];
			// printf("p = %p\n", p);
			while (*p != 'R')
				;
			// puts("check 66");
			for (int i = 0; i < count; i++)
				re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
		}
		/////////////////////////////////////////////////////////////////
		if (intra_rank == 0)
		{
			// root
			for (int c = 1; c < (intra_procn >> 2); ++c)
			{
				volatile char *p = SHM_bufs_Recv[(c << 2)];
				while (*p != 'R')
					;
				// puts("check 79");
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
			}

			// printf("%d %f \n", inter_rank, re[0]);
			// for(int c = 1;c<(intra_procn>>2);++c)
			// {
			// 	volatile char * p = SHM_bufs_Recv[(c<<2)];
			// 	for(int i = 0;i<count;i++)
			// 		((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			// 	*p = 'B';
			// }
		}
		else
		{
			// block leader
			//向上规约
			volatile char *p = SHM_bufs_Send[0];
			while (*p != 'S')
				;
			for (int i = 0; i < count; i++)
				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			*p = 'R';
			// //等待root leader结果
			// while(*p != 'B');
			// for (int i = 0; i < count; i++)
			// 	re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
			// *p = 'S';
		}
		////////////////////////////////////////////////////////////////////////////////
		//将结果广播回其它进程
		// for(int c = 1;c<=3;++c)
		// {
		// 	//对每个孩子
		// 	volatile char *p = SHM_bufs_Recv[intra_rank + c];
		// 	for(int i = 0;i<count;i++)
		// 		((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
		// 	*p = 'B';
		// }
	}
	else
	{
		volatile char *p = SHM_bufs_Send[block_leader];
		while (*p != 'S')
			;
		// puts("check 119");
		for (int i = 0; i < count; i++)
			((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
		*p = 'R';
		//等待广播结果
		// while(*p != 'B');
		// for (int i = 0; i < count; i++)
		// 	re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
		// *p = 'S';
	}
	// __sync_synchronize();
}
void one_double_allreduce_intra_BCAST_from_0(void *sendbuf, void *recvbuf, int count)
{
	// __sync_synchronize();
	// if(inter_rank == 0 && intra_rank == 0)
	// 	puts("check bcast");
	// double * re = sendbuf;
	// //puts("check");
	// if (_GLEXCOLL.intra_rank == 0)
	// {
	// 	for (int target = 0; target < ppn - 1; ++target)
	// 	{
	// 		volatile char *p = (_GLEXCOLL.SHM_buf + (target << 6));
	// 		for (int i = 0; i < count; i++)
	// 			((volatile SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
	// 		*p = 'B';
	// 	}
	// }
	// else
	// {
	// 	volatile char *p = _GLEXCOLL.SHM_buf;
	// 	while (*p != 'B')
	// 		;
	// 	for (int i = 0; i < count; i++)
	// 		((double *)recvbuf)[i] = ((volatile SHM_FLIT *)(p))->payload.T.data_double[i];
	// 	*p = 'S';
	// }
	int block_leader = (intra_rank >> 2) << 2;
	double *re = (double *)recvbuf;
	if (block_leader == intra_rank)
	{
		/////////////////////////////////////////////////////////////////
		if (intra_rank == 0)
		{
			// root
			//  //printf("%f \n",re[0]);
			// printf("%d \n",_GLEXCOLL.global_rank);
			for (int c = 1; c < (intra_procn >> 2); ++c)
			{
				volatile char *p = SHM_bufs_Recv[(c << 2)];
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p)->payload.T.data_double[i] = ((double *)sendbuf)[i];
				*p = 'B';
			}
			for (int i = 0; i < count; i++)
				((double *)re)[i] = ((double *)sendbuf)[i];
		}
		else
		{
			// block leader
			volatile char *p = SHM_bufs_Send[0];
			//等待root leader结果
			while (*p != 'B')
				;
			for (int i = 0; i < count; i++)
				re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
			*p = 'S';
		}
		////////////////////////////////////////////////////////////////////////////////
		//将结果广播回其它进程
		for (int c = 1; c <= 3; ++c)
		{
			//对每个孩子
			volatile char *p = SHM_bufs_Recv[intra_rank + c];
			for (int i = 0; i < count; i++)
				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			*p = 'B';
		}
	}
	else
	{
		volatile char *p = SHM_bufs_Send[block_leader];
		//等待广播结果
		while (*p != 'B')
			;
		for (int i = 0; i < count; i++)
			re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
		*p = 'S';
	}
}

void one_double_allreduce_intra_REDUCE_BCAST_QUEUE_sim(void *sendbuf, void *recvbuf, int count)
{
	double *re = recvbuf;
	if (intra_rank == 0)
	{ // root
		for (int i = 0; i < count; ++i)
			re[i] = ((double *)sendbuf)[i];
		// sleep(1);
		for (int i = 1; i < intra_procn; ++i)
		{
			volatile char *p = SHM_bufs_Recv[i];
			volatile char *p1 = SHM_bufs_Send[i];
			while (*p != 'R')
				;
			__sync_synchronize();
			// printf("check *p =%c\n",*p);
			// fflush(stdout);
			for (int i = 0; i < count; i++)
				re[i] += ((SHM_FLIT *)p1)->payload.T.data_double[i];
			__sync_synchronize();
			*p = 'S';
		}
		// puts("check 128");
		for (int target = 1; target < intra_procn; ++target)
		{
			volatile char *p = SHM_bufs_Recv[target];
			volatile char *p1 = SHM_bufs_Send[target];
			while (*p != 'S')
				;
			__sync_synchronize();
			for (int i = 0; i < count; i++)
				((volatile SHM_FLIT *)p1)->payload.T.data_double[i] = re[i];
			__sync_synchronize();
			*p = 'B';
		}
	}
	else
	{
		volatile char *p = SHM_bufs_Send[0];
		volatile char *p1 = SHM_bufs_Recv[0];
		while (*p != 'S')
			;
		__sync_synchronize();
		((volatile SHM_FLIT *)p1)->payload.size = count;
		for (int i = 0; i < count; i++)
			((volatile SHM_FLIT *)p1)->payload.T.data_double[i] = ((double *)sendbuf)[i];
		__sync_synchronize();
		*p = 'R';

		while (*p != 'B')
			;
		__sync_synchronize();
		// puts("check 145");
		for (int i = 0; i < count; i++)
			((double *)recvbuf)[i] = ((volatile SHM_FLIT *)(p1))->payload.T.data_double[i];
		__sync_synchronize();
		*p = 'S';
	}
}
void one_double_allreduce_intra_REDUCE_BCAST(void *sendbuf, void *recvbuf, int count)
{

	double *re = recvbuf;
	if (intra_rank == 0)
	{ // root
		for (int i = 0; i < count; ++i)
			re[i] = ((double *)sendbuf)[i];
		// sleep(1);
		for (int i = 0; i < intra_procn - 1; ++i)
		{
			volatile char *p = SHM_bufs_Recv[i + 1];
			while (*p != 'R')
				;
			// printf("check *p =%c\n",*p);
			// fflush(stdout);
			for (int i = 0; i < count; i++)
				re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
		}
		// puts("check 128");
		for (int target = 0; target < intra_procn - 1; ++target)
		{
			volatile char *p = SHM_bufs_Recv[target + 1];
			// while(*p != 'S');
			for (int i = 0; i < count; i++)
				((volatile SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			*p = 'B';
		}
	}
	else
	{
		volatile char *p = SHM_bufs_Send[0];
		while (*p != 'S')
			;
		((volatile SHM_FLIT *)p)->payload.size = count;
		for (int i = 0; i < count; i++)
			((volatile SHM_FLIT *)p)->payload.T.data_double[i] = ((double *)sendbuf)[i];
		*p = 'R';
		while (*p != 'B')
			;
		// puts("check 145");
		for (int i = 0; i < count; i++)
			((double *)recvbuf)[i] = ((volatile SHM_FLIT *)(p))->payload.T.data_double[i];
		*p = 'S';
	}

	// if(intra_rank == 0)
	// 	puts("check");
}
void SHM_Barrier(volatile int *p)
{
	//
	// while(*p !=intra_procn);
	// __sync_sub_and_fetch(p,1);
	// while(*p != 0);
	if (intra_rank != 0)
	{
		__sync_add_and_fetch(p, 1);
		// printf("*p = %d intra_procn - 1 =%d\n",*p,intra_procn - 1);
		while (*p != 0)
			;
	}
	else
	{
		while (__sync_bool_compare_and_swap(p, intra_procn - 1, 0))
		{
		}
		// printf("*p = %d\n",*p);
		// puts("release");
	}
}

void one_double_allreduce_intra_recursive_doubling(void *sendbuf, void *recvbuf, int count)
{
	double *re = (double *)recvbuf;
	for (int i = 0; i < count; ++i)
		re[i] = ((double *)sendbuf)[i];
	int step = 1;
	while (step < intra_procn)
	{
		int t = intra_rank / step;
		int target;
		if (t % 2 == 0)
		{
			target = intra_rank + step;
			// left进程
			// printf("step = %d left %d -> right = %d\n",step,intra_rank,target);
		}
		else
		{
			// right进程
			target = intra_rank - step;
			// printf("step = %d right %d->left = %d\n",step,intra_rank,target);
		}
		//先把数据写到target
		volatile char *p = SHM_bufs_Send[target]; //_GLEXCOLL.SHM_bufs[target];
		// p = p+(intra_rank<<6);
		while (*p != 'S')
			;
		for (int i = 0; i < count; i++)
			((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
		*p = 'R';
		//再等待target把数据传输来

		// p = _GLEXCOLL.SHM_bufs[intra_rank];
		// p = p+(target<<6);
		p = SHM_bufs_Recv[target];
		while (*p != 'R')
			;
		for (int i = 0; i < count; i++)
			re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
		*p = 'S';

		step = step << 1;
	}
	// 	printf("%f\n",((double *)recvbuf)[0]);
}
void one_double_allreduce_intra_l2_cache_performance_aware(void *sendbuf, void *recvbuf, int count)
{ //第一步：进行同L2cache的进程规约
	static int latency_vec[] = {20, 16, 24, 28, 8, 12, 4};

	int block_leader = (intra_rank >> 2) << 2;
	double *re = (double *)recvbuf;
	for (int i = 0; i < count; i++)
		re[i] = ((double *)sendbuf)[i];
	if (block_leader == intra_rank)
	{
		//第一步规约到block_leader
		for (int c = 1; c <= 3; ++c)
		{
			//对每个孩子
			volatile char *p = SHM_bufs_Recv[intra_rank + c];
			while (*p != 'R')
				;
			// puts("check 283");
			for (int i = 0; i < count; i++)
				re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
		}
		// printf("%d %f\n",intra_rank,re[0]);
		/////////////////////////////////////////////////////////////////
		if (intra_rank == 16)
		{
			// root
			// for (int c = 1; c < (intra_procn >> 2); ++c)
			// for (int c = 6; c >=0; --c)
			{
				volatile char *p1 = SHM_bufs_Recv[0];
				volatile char *p2 = SHM_bufs_Recv[4];
				volatile char *p3 = SHM_bufs_Recv[8];
				volatile char *p4 = SHM_bufs_Recv[12];
				volatile char *p5 = SHM_bufs_Recv[20];
				volatile char *p6 = SHM_bufs_Recv[24];
				volatile char *p7 = SHM_bufs_Recv[28];
				while (*p1 != 'R')
					;
				while (*p2 != 'R')
					;
				while (*p3 != 'R')
					;
				while (*p4 != 'R')
					;
				while (*p5 != 'R')
					;
				while (*p6 != 'R')
					;
				while (*p7 != 'R')
					;
				// while (*p7 != 'R');
				// while (*p6 != 'R');
				// while (*p5 != 'R');
				// while (*p4 != 'R');
				// while (*p3 != 'R');
				// while (*p2 != 'R');
				// while (*p1 != 'R');
				// puts("check 296");
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p1)->payload.T.data_double[i];
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p2)->payload.T.data_double[i];
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p3)->payload.T.data_double[i];
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p4)->payload.T.data_double[i];
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p5)->payload.T.data_double[i];
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p6)->payload.T.data_double[i];
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p7)->payload.T.data_double[i];
				//开始广播
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p7)->payload.T.data_double[i] = re[i];
				*p7 = 'B';
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p6)->payload.T.data_double[i] = re[i];
				*p6 = 'B';
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p5)->payload.T.data_double[i] = re[i];
				*p5 = 'B';
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p4)->payload.T.data_double[i] = re[i];
				*p4 = 'B';
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p3)->payload.T.data_double[i] = re[i];
				*p3 = 'B';
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p2)->payload.T.data_double[i] = re[i];
				*p2 = 'B';
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p1)->payload.T.data_double[i] = re[i];
				*p1 = 'B';
			}
		}
		else
		{
			// block leader
			//向上规约
			volatile char *p = SHM_bufs_Send[16];
			while (*p != 'S')
				;
			for (int i = 0; i < count; i++)
				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			*p = 'R';
			//等待root leader结果
			while (*p != 'B')
				;
			for (int i = 0; i < count; i++)
				re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
			*p = 'S';
		}
		////////////////////////////////////////////////////////////////////////////////
		//将结果广播回其它进程
		for (int c = 1; c <= 3; ++c)
		{
			//对每个孩子
			volatile char *p = SHM_bufs_Recv[intra_rank + c];
			for (int i = 0; i < count; i++)
				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			*p = 'B';
		}
	}
	else
	{
		volatile char *p = SHM_bufs_Send[block_leader];
		while (*p != 'S')
			;
		for (int i = 0; i < count; i++)
			((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
		*p = 'R';
		//等待广播结果
		while (*p != 'B')
			;
		for (int i = 0; i < count; i++)
			re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
		*p = 'S';
	}
}
void one_double_allreduce_intra_l2_cache_aware(void *sendbuf, void *recvbuf, int count)
{
	//第一步：进行同L2cache的进程规约
	int block_leader = (intra_rank >> 2) << 2;
	double *re = (double *)recvbuf;
	for (int i = 0; i < count; i++)
		re[i] = ((double *)sendbuf)[i];
	if (block_leader == intra_rank)
	{
		//第一步规约到block_leader
		for (int c = 1; c <= 3; ++c)
		{
			//对每个孩子
			volatile char *p = SHM_bufs_Recv[intra_rank + c];
			while (*p != 'R')
				;
			// puts("check 283");
			for (int i = 0; i < count; i++)
				re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
		}
		// printf("%d %f\n",intra_rank,re[0]);
		/////////////////////////////////////////////////////////////////
		if (intra_rank == 0)
		{
			// root
			for (int c = 1; c < (intra_procn >> 2); ++c)
			{
				volatile char *p = SHM_bufs_Recv[(c << 2)];
				while (*p != 'R')
					;
				// puts("check 296");
				for (int i = 0; i < count; i++)
					re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
			}
			// printf("%f \n",re[0]);
			for (int c = 1; c < (intra_procn >> 2); ++c)
			{
				volatile char *p = SHM_bufs_Recv[(c << 2)];
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
				*p = 'B';
			}
		}
		else
		{
			// block leader
			//向上规约
			volatile char *p = SHM_bufs_Send[0];
			while (*p != 'S')
				;
			for (int i = 0; i < count; i++)
				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			*p = 'R';
			//等待root leader结果
			while (*p != 'B')
				;
			for (int i = 0; i < count; i++)
				re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
			*p = 'S';
		}
		////////////////////////////////////////////////////////////////////////////////
		//将结果广播回其它进程
		for (int c = 1; c <= 3; ++c)
		{
			//对每个孩子
			volatile char *p = SHM_bufs_Recv[intra_rank + c];
			for (int i = 0; i < count; i++)
				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
			*p = 'B';
		}
	}
	else
	{
		volatile char *p = SHM_bufs_Send[block_leader];
		while (*p != 'S')
			;
		for (int i = 0; i < count; i++)
			((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
		*p = 'R';
		//等待广播结果
		while (*p != 'B')
			;
		for (int i = 0; i < count; i++)
			re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
		*p = 'S';
	}

	// if(intra_rank == 31)
	// 	printf("%f\n",re[0]/32.0);

	// if(intra_rank == 31)
	// 	printf("%f\n",re[0]/32.0);
	// if(block_leader == intra_rank)
	// {
	// 	//block leader
	// 	//printf("block leader %d \n",intra_rank);
	// 	for(int i = 0;i<count;i++) re[i] = ((double *)sendbuf)[i];

	// 	volatile char *p = (char *)_GLEXCOLL.SHM_bufs[intra_rank];
	// 	p+= (intra_rank+1)<<6;
	// 	//对每一个孩子
	// 	for(int c = 0;c<3;c++)
	// 	{
	// 		while(*p != 'R');
	// 		for(int i = 0;i<count;i++)
	// 			re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
	// 		*p = 'S';
	// 		p+=64;
	// 	}

	// 	if(intra_rank == 0)
	// 	{
	// 		//root
	// 		p =  (char *)_GLEXCOLL.SHM_bufs[0];
	// 		for(int i = 4;i<intra_procn - 1;i+=4)
	// 		{
	// 			volatile char *p1 = p + (i<<6);
	// 			while (*p1 != 'R');
	// 			for(int i = 0;i<count;i++)
	// 				re[i] += ((SHM_FLIT *)p1)->payload.T.data_double[i];
	// 			*p1 = 'S';
	// 		}
	// 		//printf("%f \n",re[0]);

	// 		//接下来root开始广播
	// 		for(int target = 1;target<intra_procn ;++target)
	// 		{
	// 			volatile char *p = _GLEXCOLL.SHM_bufs[target];
	// 			for(int i = 0;i<count;i++)
	// 				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
	// 			*p = 'B';
	// 		}
	// 	}else{
	// 		//block leader
	// 		//将块写到root上去
	// 		p = (char *)_GLEXCOLL.SHM_bufs[0];
	// 		p+= (intra_rank<<6);
	// 		while (*p!='S');
	// 		for(int i = 0;i<count;i++)
	// 			((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
	// 		*p = 'R';
	// 		//等待root广播结果
	// 		p = (char *)_GLEXCOLL.SHM_bufs[intra_rank];
	// 		while(*p != 'B');
	// 			for(int i = 0;i<count;i++)
	// 				re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
	// 		*p = 'S';
	// 	}
	// }else{
	// 	//block leaf
	// 	//将块写到block leader上去。
	// 	volatile char * p =(volatile char *) _GLEXCOLL.SHM_bufs[block_leader];
	// 	p += (intra_rank<<6);
	// 	while(*p != 'S');
	// 	for(int i = 0;i<count;i++)
	// 		((SHM_FLIT *)p)->payload.T.data_double[i] = ((double *)sendbuf)[i];
	// 	*p = 'R';

	// 	//等待root广播结果
	// 	p = (char *)_GLEXCOLL.SHM_bufs[intra_rank];
	// 	while (*p != 'B');
	// 	for (int i = 0; i < count; i++)
	// 		re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
	// 	*p = 'S';
	// }
}
void one_double_allreduce_intra_2_NOMINAL(void *sendbuf, void *recvbuf, int count)
{
	int n = intra_procn;
	double *re = (double *)recvbuf;
	for (int i = 0; i < count; ++i)
		re[i] = ((double *)sendbuf)[i];
	int step = 1;
	//规约过程

	while (step <= n - 1)
	{
		if (intra_rank % step == 0)
		{
			int p = intra_rank / step;
			if (p % 2 == 0)
			{
				int source = intra_rank + step;
				if (source < intra_procn)
				{
					//接收进程
					// printf("step = %d intra_rank  = %d recv\n",step,intra_rank);
					// volatile char *p=(volatile char *)_GLEXCOLL.SHM_bufs[intra_rank];
					// p+=(source<<6);
					volatile char *p = SHM_bufs_Recv[source];
					while (*p != 'R')
						;
					for (int i = 0; i < count; i++)
						re[i] += ((SHM_FLIT *)p)->payload.T.data_double[i];
				}
			}
			else
			{
				//发送进程
				int target = intra_rank - step;
				// printf("step = %d intra_rank  = %d Send\n",step,intra_rank);
				//  volatile char *p=(volatile char *)_GLEXCOLL.SHM_bufs[target];
				//  p += (intra_rank<<6);
				volatile char *p = SHM_bufs_Send[target];
				while (*p != 'S')
					;
				for (int i = 0; i < count; i++)
					((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
				*p = 'R';
			}
		}
		step = (step << 1);
	}
	// printf("%d %f\n",intra_rank,re[0]);
	//接下来是广播过程。
	step = (step >> 1);
	while (step > 0)
	{
		if (intra_rank % step == 0)
		{
			int p = intra_rank / step;
			if (p % 2 == 0)
			{
				int target = intra_rank + step;
				if (target < intra_procn)
				{
					//发送进程
					// printf("step = %d intra_rank  = %d recv\n",step,intra_rank);
					// volatile char *p=(volatile char *)_GLEXCOLL.SHM_bufs[intra_rank];
					// p+=(source<<6);
					volatile char *p = SHM_bufs_Recv[target];
					for (int i = 0; i < count; i++)
						((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
					*p = 'B';
				}
			}
			else
			{
				//接收进程
				int source = intra_rank - step;
				// printf("step = %d intra_rank  = %d Send\n",step,intra_rank);
				//  volatile char *p=(volatile char *)_GLEXCOLL.SHM_bufs[target];
				//  p += (intra_rank<<6);
				volatile char *p = SHM_bufs_Send[source];
				while (*p != 'B')
					;
				for (int i = 0; i < count; i++)
					re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
				*p = 'S';
			}
		}
		step = (step >> 1);
	}

	// MPI_Barrier(Comm_intra);
	// if(intra_rank == 0){
	// 		for(int target = 1;target<intra_procn ;++target)
	// 		{
	// 			volatile char *p = _GLEXCOLL.SHM_bufs[target];
	// 			for(int i = 0;i<count;i++)
	// 				((SHM_FLIT *)p)->payload.T.data_double[i] = re[i];
	// 			*p = 'B';
	// 		}
	//  }else{
	// 	volatile char *p = (char *)_GLEXCOLL.SHM_bufs[intra_rank];
	// 	while (*p != 'B');
	// 	for (int i = 0; i < count; i++)
	// 		re[i] = ((SHM_FLIT *)p)->payload.T.data_double[i];
	// 	*p = 'S';
	//  }
	// if(intra_rank == 23)
	// 	printf("%f\n",re[0]);
}
void one_double_allreduce_intra(void *sendbuf, void *recvbuf, int count)
{
	switch (INTRA_ALLREDUCE_TYPE)
	{
	case 0:
		/* REDUCE_BCAST */
		one_double_allreduce_intra_REDUCE_BCAST(sendbuf, recvbuf, count);
		break;
	case 1:
		/*2_NOMINAL*/
		one_double_allreduce_intra_2_NOMINAL(sendbuf, recvbuf, count);
		break;
	case 2:
		/*Recurseive doubling*/
		one_double_allreduce_intra_recursive_doubling(sendbuf, recvbuf, count);
		break;
	case 3:
		/*Recurseive doubling*/
		one_double_allreduce_intra_l2_cache_aware(sendbuf, recvbuf, count);
		break;
	case 4:
		/*Intel REDUCE_BCAST reorder*/
		one_double_allreduce_intra_REDUCE_BCAST_Intel(sendbuf, recvbuf, count);
		break;
	case 5:
		/**/
		one_double_allreduce_intra_REDUCE_BCAST_QUEUE_sim(sendbuf, recvbuf, count);
		break;
	case 6:
		/**/
		one_double_allreduce_intra_l2_cache_performance_aware(sendbuf, recvbuf, count);
		break;

	default:
		break;
	}
}

//一下是针对更大消息的reduce和bcast
extern void *allreduce_intra_node_buf_cacheleaderBUF_header;
extern void *large_allreduce_bufferP;
extern int allreduce_CorePerCache;
extern volatile void *allreduce_my_state_on_parent;
extern volatile void *allreduce_cacheL_state_on_root;
extern volatile void *allreduce_cacheL_flagsstart_on_root;
extern void *large_allreduce_bufferROOT;
// struct sync_lock
// {
// 	int lock0[16];
// };
// #define allreduce_intra_node_locks_n 16
// struct sync_header
// {
// 	struct sync_lock mlocks[allreduce_intra_node_locks_n];
// };
#define compiler_barrier() __asm__ __volatile__("" \
												:  \
												:  \
												: "memory")

// void GLEX_Small_message_reduce_double_sum(double *sendbuf, double *recvbuf, int count)
// {
// 	int allreduce_intra_rank = intra_rank;
// 	int allreduce_intra_procn = intra_procn;
// 	int cache_leader = (allreduce_intra_rank / allreduce_CorePerCache) * allreduce_CorePerCache;
// 	if (allreduce_intra_rank != cache_leader)
// 	{
// 		//第一步等待缓冲区为空闲状态
// 		while (*(volatile int *)allreduce_my_state_on_parent != 0)
// 			;
// 		//开始拷贝数据到缓冲区,按照64字节对齐
// 		volatile double *p = (volatile double *)(large_allreduce_bufferP + ((((allreduce_intra_rank - cache_leader) * count * 8) >> 6) << 6) + 64);
// 		// printf("shift = %d\n", (void *)p - large_allreduce_bufferP);
// 		for (int i = 0; i < count; i++)
// 			p[i] = sendbuf[i];
// 		// memcpy(p, sendbuf, count * 8);
// 		//数据拷贝完成后设置标志为数据可用状态
// 		compiler_barrier();
// 		mfence();
// 		*(volatile int *)allreduce_my_state_on_parent = 996;
// 	}
// 	else
// 	{
// 		// memcpy(recvbuf, sendbuf, count * 8);
// 		//cache leader节点负责将数据累加到recv缓冲区中
// 		for (int c = 1; c < allreduce_CorePerCache; c++)
// 		{
// 			volatile int *waitp = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[c].lock0;
// 			while (*waitp != 996)
// 				;
// 			volatile double *recvp = (double *)(large_allreduce_bufferP + (((c * count * 8) >> 6) << 6) + 64);
// 			// printf("shift = %d\n", (void *)recvp - large_allreduce_bufferP);
// 			compiler_barrier();
// 			mfence();
// 			for (int i = 0; i < count; i++)
// 			{
// 				if (c == 1)
// 					((double *)large_allreduce_bufferP)[i] = sendbuf[i] + recvp[i];
// 				else
// 					((double *)large_allreduce_bufferP)[i] += recvp[i];
// 			}
// 			*waitp = 0;
// 		}
// 		if (*(volatile double *)large_allreduce_bufferP < 3.99 || *(volatile double *)large_allreduce_bufferP > 4.001)
// 		{
// 			printf("%d %d error reduce up %f\n", intra_rank, inter_rank, *(volatile double *)large_allreduce_bufferP);
// 		}
// 		// printf("large_allreduce_bufferP %f\n", *(double *)large_allreduce_bufferP);
// 		//第二步继续向上规约。
// 		{
// 			if (allreduce_intra_rank != 0)
// 			{
// 				// cache leaders
// 				int c = allreduce_intra_rank / allreduce_CorePerCache;
// 				double *buffer_on_root = (double *)(allreduce_cacheL_flagsstart_on_root + sizeof(struct sync_header)) + ((((c + allreduce_CorePerCache) * count * 8) >> 6) << 6) + 64;
// 				while (*(volatile int *)allreduce_cacheL_state_on_root != 0)
// 					;
// 				//将数据拷贝进去
// 				for (int i = 0; i < count; i++)
// 				{
// 					buffer_on_root[i] = ((double *)large_allreduce_bufferP)[i];
// 				}
// 				*(volatile int *)allreduce_cacheL_state_on_root = 997;
// 			}
// 			else
// 			{
// 				int childn = allreduce_intra_procn / allreduce_CorePerCache;
// 				for (int c = 1; c < childn; c++)
// 				{
// 					volatile int *waitp = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[c + allreduce_CorePerCache].lock0;
// 					volatile double *datap = (volatile double *)large_allreduce_bufferP + ((((c + allreduce_CorePerCache) * count * 8) >> 6) << 6) + 64;
// 					//  double *datap = ( double *)(allreduce_intra_node_buf_cacheleaderBUF_header + sizeof(allreduce_Header)) + ((((c + allreduce_CorePerCache) * count * 8) >> 6) << 6) + 64;
// 					// if(allreduce_intra_rank == 0)
// 					//         int test = (unsigned long long)large_allreduce_bufferP - (unsigned long long)datap;
// 					// printf("%p , %p,%p\n",allreduce_intra_node_buf_cacheleaderBUF_header,large_allreduce_bufferP,datap);
// 					while (*waitp != 997)
// 					{
// 						/* code */
// 					}
// 					// lfence();
// 					// printf("%f\n", *(volatile double *)datap);
// 					// printf("%f+%f\n",*(volatile double *)large_allreduce_bufferP, *(volatile double *)datap);

// 					for (int i = 0; i < count; i++)
// 					{
// 						if (c == 1)
// 						{
// 							recvbuf[i] = ((volatile double *)large_allreduce_bufferP)[i] + datap[i];
// 						}
// 						else
// 						{
// 							recvbuf[i] += datap[i];
// 						}
// 					}
// 					// sfence();
// 					*waitp = 0;
// 				}
// 				// printf("%f\n", *(volatile double *)recvbuf);
// 				//root 进程
// 			}
// 		}
// 	}
// }

#define compile_barrier() __asm__ __volatile__("" \
											   :  \
											   :  \
											   : "memory");

#pragma GCC push_options
#pragma GCC optimze("O3")
void GLEX_Small_message_reduce_double_sum(double *sendbuf, double *recvbuf, int count)
{
	int allreduce_intra_rank = intra_rank;
	int allreduce_intra_procn = intra_procn;
	int cache_leader = (allreduce_intra_rank / allreduce_CorePerCache) * allreduce_CorePerCache;
	if (allreduce_intra_rank != cache_leader)
	{
		//第一步等待缓冲区为空闲状态
		while (*(volatile int *)allreduce_my_state_on_parent != 0)
			;
		compile_barrier();
		//开始拷贝数据到缓冲区,按照64字节对齐
		volatile double *p = (double *)(large_allreduce_bufferP + ((((allreduce_intra_rank - cache_leader) * count * 8) >> 6) << 6) + 64);
// printf("shift = %d\n", (void *)p - large_allreduce_bufferP);
#pragma omp simd
		for (int i = 0; i < count; i++)
			p[i] = sendbuf[i];
		// memcpy(p, sendbuf, count * 8);
		//数据拷贝完成后设置标志为数据可用状态
		__sync_synchronize();
		compile_barrier();
		*(volatile int *)allreduce_my_state_on_parent = 996;
		//等待父节点接收消息完毕
		while (*(volatile int *)allreduce_my_state_on_parent != 0)
			;
	}
	else
	{
		// MPI_Barrier(Comm_intra);

		// memcpy(recvbuf, sendbuf, count * 8);
		// cache leader节点负责将数据累加到recv缓冲区中
		for (int c = 1; c < allreduce_CorePerCache; c++)
		{
			volatile int *waitp = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[c].lock0;
			compile_barrier();
			while (*waitp != 996)
				;
			compile_barrier();
			__sync_synchronize();
			double *recvp = (double *)(large_allreduce_bufferP + (((c * count * 8) >> 6) << 6) + 64);
			// if (*(volatile double *)recvp < 0.99 || *(volatile double *)recvp > 1.001)
			// {
			// 	printf("%d error recvp  %f\n", allreduce_intra_rank, *(volatile double *)recvp);
			// }
			// puts("check 1150");
			// printf("shift = %d\n", (void *)recvp - large_allreduce_bufferP);
			// if (intra_rank == 0)
			if (c == 1)
			{
#pragma omp simd
				for (int i = 0; i < count; i++)
					((volatile double *)large_allreduce_bufferP)[i] = sendbuf[i] + recvp[i];
			}
			else
			{
#pragma omp simd
				for (int i = 0; i < count; i++)
				{
					((volatile double *)large_allreduce_bufferP)[i] += recvp[i];
				}
			}

			compile_barrier();
			*waitp = 0;
		}

		// printf("rank=%d reduce re: %lf\n", global_rank, ((volatile double *)large_allreduce_bufferP)[0]);
		// if (*(volatile double *)large_allreduce_bufferP < 3.99 || *(volatile double *)large_allreduce_bufferP > 4.001)
		// {
		// 	printf("%d error reduce up %f\n", allreduce_intra_rank, *(volatile double *)large_allreduce_bufferP);
		// }
		// printf("large_allreduce_bufferP %f\n", *(double *)large_allreduce_bufferP);
		//第二步继续向上规约。
		if (allreduce_intra_procn > allreduce_CorePerCache)
		{
			if (allreduce_intra_rank != 0)
			{
				// cache leaders
				int c = allreduce_intra_rank / allreduce_CorePerCache;
				volatile double *buffer_on_root = (volatile double *)(allreduce_cacheL_flagsstart_on_root + sizeof(struct sync_header));
				compile_barrier();
				while (*(volatile int *)allreduce_cacheL_state_on_root != 0)
					;
				compile_barrier();
				//将数据拷贝进去
				for (int i = 0; i < count; i++)
				{
					buffer_on_root[i] = ((double *)large_allreduce_bufferP)[i];
				}
				compile_barrier();
				__sync_synchronize();

				*(volatile int *)allreduce_cacheL_state_on_root = 997;
				while (*(volatile int *)allreduce_cacheL_state_on_root != 0)
					;
			}
			else
			{
				int childn = allreduce_intra_procn / allreduce_CorePerCache;
				for (int c = 1; c < childn; c++)
				{
					volatile int *waitp = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[c + allreduce_CorePerCache].lock0;
					volatile double *datap = (volatile double *)(large_allreduce_bufferP);
					compile_barrier();
					//  double *datap = ( double *)(allreduce_intra_node_buf_cacheleaderBUF_header + sizeof(allreduce_Header)) + ((((c + allreduce_CorePerCache) * count * 8) >> 6) << 6) + 64;
					// if(allreduce_intra_rank == 0)
					//         int test = (unsigned long long)large_allreduce_bufferP - (unsigned long long)datap;
					// printf("%p , %p,%p\n",allreduce_intra_node_buf_cacheleaderBUF_header,large_allreduce_bufferP,datap);
					while (*waitp != 997)
					{
						/* code */
					}
					compile_barrier();
					// lfence();
					// printf("%f\n", *(volatile double *)datap);
					// printf("%f+%f\n",*(volatile double *)large_allreduce_bufferP, *(volatile double *)datap);
					__sync_synchronize();

					for (int i = 0; i < count; i++)
					{
						if (c == 1)
						{
							recvbuf[i] = ((volatile double *)large_allreduce_bufferP)[i] + datap[i];
						}
						else
						{
							recvbuf[i] += datap[i];
						}
						// if (c == childn - 1)
						// 	printf("i = %d %d check recvb %f %f\n", i, global_rank, recvbuf[i], ((volatile double *)datap)[i]);
					}
					// __sync_synchronize();

					compile_barrier();
					// sfence();
					*waitp = 0;
				}
				// printf("%f\n", *(volatile double *)recvbuf);
				// root 进程
			}
		}
		else
		{
			memcpy(recvbuf, large_allreduce_bufferP, count * sizeof(double));
		}
	}
}

void GLEX_Small_message_bcast_double(double *sendbuf, double *recvbuf, int count)
{
	//规约和广播不能用同一个内存体系，否则会出现错误。
	int shift = (((allreduce_CorePerCache * count * 8) >> 6) << 6) + 64;
	int allreduce_intra_rank = intra_rank;
	int allreduce_intra_procn = intra_procn;
	int cache_leader = (allreduce_intra_rank / allreduce_CorePerCache) * allreduce_CorePerCache;
	if (allreduce_intra_rank == 0)
	{
		// root进程负责广播消息
		// root进程首先等待缓存空间为空闲状态
		volatile double *datap = (volatile double *)(large_allreduce_bufferP + ((((allreduce_CorePerCache)*count * 8) >> 6) << 6) + 64);
		// while (*waitp != 0)
		// 	;
		// *waitp = 999;
		if (sendbuf != recvbuf)
		{
			for (int i = 0; i < count; i++)
			{
				recvbuf[i] = datap[i] = sendbuf[i];
			}
		}
		else
		{
			for (int i = 0; i < count; i++)
			{
				datap[i] = sendbuf[i];
			}
		}
		// printf("result sent to child = %f\n", *datap);
		__sync_synchronize();
		// sleep(1);
		int cache_leaderN = (allreduce_intra_procn / allreduce_CorePerCache);
		for (int c = 1; c < cache_leaderN; c++)
		{
			volatile int *flagp = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[allreduce_CorePerCache + c].lock0;
			while (*flagp != 0)
				;
			__sync_synchronize();
			*flagp = 1000;
		}
		for (int c = 1; c < cache_leaderN; c++)
		{
			volatile int *flagp = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[allreduce_CorePerCache + c].lock0;
			while (*flagp != 1001)
				;
			__sync_synchronize();
			*flagp = 0;
		}

		// *waitp = 0;
	}
	else
	{
		if (allreduce_intra_rank % allreduce_CorePerCache == 0)
		{

			// cache leader节点需要等待父节点的消息
			volatile int *flagw = (volatile int *)allreduce_cacheL_state_on_root;
			while (*flagw != 1000)
				;
			__sync_synchronize();
			volatile double *buffer_on_parent = (volatile double *)(large_allreduce_bufferROOT + ((((allreduce_CorePerCache)*count * 8) >> 6) << 6) + 64);
			// volatile int *flagTochild = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[0];
			volatile double *buffer_to_child = (volatile double *)(large_allreduce_bufferP + ((((allreduce_CorePerCache)*count * 8) >> 6) << 6) + 64);

			for (int i = 0; i < count; i++)
			{
				recvbuf[i] = buffer_to_child[i] = buffer_on_parent[i];
			}
			__sync_synchronize();
			*flagw = 1001;
		}
	}
	// if(allreduce_intra_rank % allreduce_CorePerCache == 0)
	//     printf("intra_rank = %d result recv from root = %f\n",allreduce_intra_rank, *(double *)recvbuf);
	// Cache Leader节点向剩余进程广播数据
	// if (0)
	if (allreduce_intra_rank % allreduce_CorePerCache == 0)
	{
		// for (int i = 0; i < count; i++)
		// {
		// 	recvbuf[i] = ((volatile double *)large_allreduce_bufferP)[i];
		// }
		// printf("intra_rank = %d result recv from root = %f\n",allreduce_intra_rank, *(double *)recvbuf);
		// cache leader 节点
		for (int c = 1; c < allreduce_CorePerCache; c++)
		{
			volatile int *flagTochild = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[c].lock0;
			while (*flagTochild != 0)
				;
			__sync_synchronize();
			*flagTochild = 98;
		}
	}
	else
	{
		// child
		int shif = allreduce_intra_rank - cache_leader;
		volatile int *flag_on_parent = ((struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header)->mlocks[shif].lock0;
		while (*flag_on_parent != 98)
			;
		__sync_synchronize();
		volatile double *buffer_on_parent = (volatile double *)(large_allreduce_bufferP + ((((allreduce_CorePerCache)*count * 8) >> 6) << 6) + 64);
		for (int i = 0; i < count; i++)
		{
			recvbuf[i] = ((double *)buffer_on_parent)[i];
		}
		__sync_synchronize();
		*flag_on_parent = 0;
	}
	// printf("%d result recv frome parent = %f\n", global_rank, *recvbuf);
	// fflush(stdout);
	// MPI_Barrier(Comm_intra);
	// exit(0);
}
static int allreduce_flag_count = 0;
void medium_reduce(double *sendbuf, double *recvbuf, int count)
{
	//第一步放置数据到自己所在共享内存区域
	if (count * 8 > allreduce_buffer_size_single)
	{
		puts("消息过大，建议转为大消息规约模式");
	}
	// if (intra_rank % allreduce_CorePerCache == 0)
	//等待缓冲区状态为0
	//0->1
	while (*(allreduce_shm_flags->mlocks[intra_rank].lock0) != allreduce_flag_count)
		;
	{
		//叶子节点将数据拷贝到对应缓冲区中
		memcpy(allreduce_shm_buffer_starts[intra_rank], sendbuf, count * sizeof(double));
		__sync_synchronize();
		allreduce_shm_flags->mlocks[intra_rank].lock0[0] = 1;
	}
	// MPI_Barrier(Comm_intra);
	//将消息规约到cache leader
	//1->2
	if (intra_rank % allreduce_CorePerCache == 0)
	{
		for (int c = 1; c < allreduce_CorePerCache; c++)
		{
			while ((allreduce_shm_flags->mlocks[intra_rank + c].lock0[0]) != 1)
				;
			// usleep(5);
			// printf("allreduce_shm_flags[%d] = %d \n", intra_rank + c, *(allreduce_shm_flags->mlocks[intra_rank + c].lock0));
			__sync_synchronize();
			double *bufp = (double *)(allreduce_shm_buffer_starts[intra_rank]);
			double *sp = (double *)(allreduce_shm_buffer_starts[intra_rank + c]);
			for (int i = 0; i < count; i++)
			{
				bufp[i] += sp[i];
			}
		}
		__sync_synchronize();
		*(allreduce_shm_flags->mlocks[intra_rank].lock0) = 2;
	}
	// else
	// {
	//     int leader = (intra_rank / allreduce_CorePerCache);
	//     leader *= allreduce_CorePerCache;
	//     while (*(allreduce_shm_flags->mlocks[leader].lock0) = !2)
	//         ;
	//     __sync_synchronize();
	//     *(allreduce_shm_flags->mlocks[intra_rank].lock0) = 2;
	// }
	//将消息规约到root
	if (intra_rank % allreduce_CorePerCache == 0)
	{
		if (intra_rank == 0)
		{
			int childn = intra_procn / allreduce_CorePerCache;
			for (int c = 1; c < childn; c++)
			{
				while (*(allreduce_shm_flags->mlocks[allreduce_CorePerCache * c].lock0) != 2)
					;
				__sync_synchronize();
				double *bufp = (double *)(allreduce_shm_buffer_starts[intra_rank]);
				double *sp = (double *)(allreduce_shm_buffer_starts[allreduce_CorePerCache * c]);
				for (int i = 0; i < count; i++)
				{
					bufp[i] += sp[i];
				}
				// if (intra_rank == 0)
				// {
				//     printf("%f \n", bufp[0]);
				// }
			}
			memcpy(recvbuf, allreduce_shm_buffer_starts[intra_rank], count * sizeof(double));
		}
		// else
		// {
		//     while (*(allreduce_shm_flags->mlocks[0].lock0) != 3)
		//         ;
		// }
	}
	// __sync_synchronize();
	// *(allreduce_shm_flags->mlocks[intra_rank].lock0) = 3;
	// = 3;
}

void medium_bcast(double *sendbuf, double *recvbuf, int count)
{
	if (intra_rank % allreduce_CorePerCache == 0)
	{
		if (intra_rank == 0)
		{
			//叶子节点将数据拷贝到对应缓冲区中
			memcpy(allreduce_shm_buffer_starts[intra_rank], sendbuf, count * sizeof(double));
			__sync_synchronize();
			*(allreduce_shm_flags->mlocks[intra_rank].lock0) = 3;
			// __sync_synchronize();
			int childn = intra_procn / allreduce_CorePerCache;
			for (int c = 1; c < childn; c++)
			{
				while (*(allreduce_shm_flags->mlocks[c * allreduce_CorePerCache].lock0) != 2)
					;
				__sync_synchronize();
				memcpy(allreduce_shm_buffer_starts[c * allreduce_CorePerCache], sendbuf, count * sizeof(double));
				__sync_synchronize();
				*(allreduce_shm_flags->mlocks[c * allreduce_CorePerCache].lock0) = 3;
			}
		}
		// else
		// {
		//     while (*(allreduce_shm_flags->mlocks[0].lock0) !=  3)
		//         ;
		//     __sync_synchronize();
		//     memcpy(allreduce_shm_buffer_starts[intra_rank], allreduce_shm_buffer_starts[0], count * sizeof(double));
		//     __sync_synchronize();
		//     *(allreduce_shm_flags->mlocks[intra_rank].lock0) =  3;
		// }
		// printf("%d get %f\n", intra_rank, ((double *)allreduce_shm_buffer_starts[intra_rank])[count - 1]);
	}

	if (intra_rank % allreduce_CorePerCache == 0)
	{
		while ((*(allreduce_shm_flags->mlocks[intra_rank].lock0) != 3))
			;
		memcpy(recvbuf, allreduce_shm_buffer_starts[intra_rank], count * sizeof(double));
		for (int c = 1; c < allreduce_CorePerCache; c++)
		{
			while ((*(allreduce_shm_flags->mlocks[intra_rank + c].lock0) != 2)) //&& (*(allreduce_shm_flags->mlocks[intra_rank + c].lock0) != 1))
				;
			//     memcpy(recvbuf, allreduce_shm_buffer_starts[leader], count * sizeof(double));
		}
		*(allreduce_shm_flags->mlocks[intra_rank].lock0) = 0;
	}
	else
	{
		int leader = (intra_rank / allreduce_CorePerCache);
		leader *= allreduce_CorePerCache;
		while (*(allreduce_shm_flags->mlocks[leader].lock0) != 3)
			;
		__sync_synchronize();
		memcpy(recvbuf, allreduce_shm_buffer_starts[leader], count * sizeof(double));
		__sync_synchronize();
		*(allreduce_shm_flags->mlocks[intra_rank].lock0) = 2;

		while ((*(allreduce_shm_flags->mlocks[leader].lock0) != 0) && *(allreduce_shm_flags->mlocks[leader].lock0) != 1)
			;
		__sync_synchronize();
		*(allreduce_shm_flags->mlocks[intra_rank].lock0) = 0;
	}
	// = 4;
	// printf("%d get %f\n", intra_rank, recvbuf[count - 1]);
	// fflush(stdout);
}

#pragma GCC pop_options