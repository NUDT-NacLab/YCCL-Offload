#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <endian.h>
#include <sched.h>
#include <errno.h>
#include <mpi.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "glex.h"
#include "glexcoll.h"
#include "BasicDataStructure.h"
#include "glexallreduce.h"

struct _GLEXCOLL_Data _GLEXCOLL;
MPI_Comm Comm_intra;
int intra_procn;
int intra_rank;
MPI_Comm Comm_inter;
int inter_procn;
int inter_rank;
int global_rank;
int global_procn;

int Childn_K;

int iphmh_num;
glex_mem_handle_t *send_mhs[GLEX_MELLOC_MAX];
glex_mem_handle_t *recv_mhs[GLEX_MELLOC_MAX];
//设置线程亲和性
void GLEX_set_affinity(int rank)
{
	int s, j;
	cpu_set_t cpuset;
	pthread_t thread;
	int speid = rank;

	thread = pthread_self();

	/* Set affinity mask to include CPUs 0 to 7 */

	CPU_ZERO(&cpuset);
	CPU_SET(speid, &cpuset);

	s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if (s != 0)
	{
		puts("pthread_setaffinity_np error");
		exit(0);
	}

	//检查当前cpu的位置分配
	/* Check the actual affinity mask assigned to the thread */
	s = pthread_getaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if (s != 0)
	{
		// handle_error_en(s, "pthread_getaffinity_np");
		puts("pthread_getaffinity_np error");
		exit(0);
	}
	//设置内存分配只能放置在当前numa内
	// char cpuallowd[65];
	// printf("numa_num_task_cpus = %d\n",numa_num_task_cpus());
	// numa_set_membind
	{
		// if(_GLEXCOLL.intra_rank == 31)
		// {
		// 	pid_t pid = getpid();
		// 	char cmd[100] = {};
		// 	sprintf(cmd,"cat /proc/%d/task/%d/status\n",(int)pid,(int)pid);
		// 	system(cmd);
		// }
	}
}

int GLEXCOLL_max_imm_RDMA_Len;
int GLEXCOLL_max_PPSDSM_Len;
int GLEXCOLL_max_RDMA_Len;

int CorePerNuma = 32;
struct glex_lib_dev
{
	struct glex_device glex_dev;
	struct glex_lib_ep *ep_list;
	struct glex_device_attr attr;
	int fd; /* handle of NIO device */
};
int *host_ids;
int GLEXCOLL_Init(int argc, char **argv)
{

	char *PPNC = getenv("GLEX_COLL_PPN");
	if (PPNC == 0)
	{
		puts("GLEXCOLL_Init PPNC=0");
		exit(1);
	}
	sscanf(PPNC, "%d", &ppn);
	// printf("PPN = %d\n", ppn);
	// puts("start init");
	_GLEXCOLL.Comm_Global = MPI_COMM_WORLD;
	glex_ret_t ret;
	glex_num_of_device(&(_GLEXCOLL.num_of_devices));

	ret = glex_open_device(0, &(_GLEXCOLL.dev));
	if (ret != GLEX_SUCCESS)
	{
		printf("_open_device() error, return: %d\n", ret);
		exit(1);
	}
	//查询设备属性
	ret = glex_query_device(_GLEXCOLL.dev, &(_GLEXCOLL.dev_attr));
	if (ret != GLEX_SUCCESS)
	{
		printf("_query_device(), return: %d\n", ret);
		exit(1);
	}
	GLEXCOLL_max_imm_RDMA_Len = _GLEXCOLL.dev_attr.max_imm_rdma_data_len;
	GLEXCOLL_max_PPSDSM_Len = 56;
	GLEXCOLL_max_RDMA_Len = _GLEXCOLL.dev_attr.max_rdma_data_len;
	// if (global_rank == 0)
	// {
	// 	printf("max_imm_rdma_data_len = %d max_rdma_data_len=%d\n", _GLEXCOLL.dev_attr.max_imm_rdma_data_len, _GLEXCOLL.dev_attr.max_rdma_data_len);
	// 	exit(0);
	// }
	// printf("max fast ep num = %d,max ep num = %d, max emulate ep num = %d\n",_GLEXCOLL.dev_attr.max_fast_ep_num,_GLEXCOLL.dev_attr.max_ep_num,_GLEXCOLL.dev_attr.max_emulate_ep_num);
	// max fast ep num = 7,max ep num = 39, max emulate ep num = 16384
	// print("mmt units = %d\n",_GLEXCOLL.dev_attr.mmt_units);
	//创建节点内通信子
	MPI_Comm_size(_GLEXCOLL.Comm_Global, &(_GLEXCOLL.global_procn));
	// printf("(_GLEXCOLL.global_procn = %d\n",_GLEXCOLL.global_procn);
	MPI_Comm_rank(_GLEXCOLL.Comm_Global, &(_GLEXCOLL.global_rank));

	// split communicator
	int namelen, n;
	MPI_Get_processor_name(host_name, &namelen);

	// char(*host_names)[MPI_MAX_PROCESSOR_NAME];
	// //printf("host_name =%s \n",host_name);

	// int bytes = (_GLEXCOLL.global_procn) * sizeof(char[MPI_MAX_PROCESSOR_NAME]);
	// host_names = (char(*)[MPI_MAX_PROCESSOR_NAME])malloc(bytes);

	// strcpy(host_names[(_GLEXCOLL.global_rank)], host_name);
	// MPI_Allgather(host_name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
	// 				host_names,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
	// 					MPI_COMM_WORLD);
	// // for (n = 0; n < (_GLEXCOLL.global_procn); n++)
	// // 	MPI_Bcast(&(host_names[n]), MPI_MAX_PROCESSOR_NAME, MPI_CHAR, n, _GLEXCOLL.Comm_Global);
	// //按照节点的顺序重排整个进程主机名字，这样相同节点的进程就在数组中紧挨着
	// qsort(host_names, (_GLEXCOLL.global_procn), sizeof(char[MPI_MAX_PROCESSOR_NAME]), stringCmp);

	// //同一个节点内的进程划分为节点内通信子
	// //找到节点内进程数量和ID
	// int color = 0;
	// for (n = 0; n < (_GLEXCOLL.global_procn); n++)
	// {
	// 	//为同一个节点内的host划分同color颜色。
	// 	//第一个节点内的进程颜色为0,第二个节点内进程颜色为1，以此类推
	// 	if (n > 0 && strcmp(host_names[n - 1], host_names[n]))
	// 		color++;
	// 	if (strcmp(host_name, host_names[n]) == 0)
	// 		break;
	// }
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(global_rank == 0)
	// {
	// 	puts("==========================181===================================");
	// 	fflush(stdout);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);

	// MPI_Barrier(MPI_COMM_WORLD);
	// if(_GLEXCOLL.global_rank == 0)
	// {
	// 	puts("==========================181===================================");
	// 	fflush(stdout);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	int color = _GLEXCOLL.global_rank / ppn;
	MPI_Comm_split(_GLEXCOLL.Comm_Global, color, (_GLEXCOLL.global_rank), &(_GLEXCOLL.Comm_intra));
	MPI_Comm_size(_GLEXCOLL.Comm_intra, &(_GLEXCOLL.intra_procn));
	MPI_Comm_rank(_GLEXCOLL.Comm_intra, &(_GLEXCOLL.intra_rank));

	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(_GLEXCOLL.global_rank == 0)
	// {
	// 	puts("==========================194===================================");
	// 	fflush(stdout);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// printf("intra procn = %d,intra rank = %d\n", _GLEXCOLL.intra_procn, _GLEXCOLL.intra_rank);
	//现在为每个节点上loca_rank相同的进程创建Comm_inter；
	MPI_Comm_split(_GLEXCOLL.Comm_Global, _GLEXCOLL.intra_rank, (_GLEXCOLL.global_rank), &(_GLEXCOLL.Comm_inter));
	MPI_Comm_size(_GLEXCOLL.Comm_inter, &(_GLEXCOLL.inter_procn));
	MPI_Comm_rank(_GLEXCOLL.Comm_inter, &(_GLEXCOLL.inter_rank));
	//根据进程ID设置主线程亲和性：
	GLEX_set_affinity(_GLEXCOLL.intra_rank);

	// MPI_Barrier(MPI_COMM_WORLD);
	// if(_GLEXCOLL.global_rank == 0)
	// {
	// 	puts("==========================213===================================");
	// 	fflush(stdout);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	if (_GLEXCOLL.intra_rank == 0)
	{
		host_ids = (int *)malloc(sizeof(int) * _GLEXCOLL.inter_procn);
		int my_hostid = -1;
		sscanf(host_name, "cn%d", &(my_hostid));
		MPI_Allgather(&my_hostid, 1, MPI_INT, host_ids, 1, MPI_INT, _GLEXCOLL.Comm_inter);
	}
	// free(host_names);
	//根据每个numa16核心为分配numa intra communicator
	color = _GLEXCOLL.intra_rank / CorePerNuma;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(global_rank == 0)
	// {
	// 	puts("199");
	// 	fflush(stdout);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// printf("grank = %d, intra_rank = %d, color=%d CorePerNuma = %d\n", (_GLEXCOLL.global_rank), _GLEXCOLL.intra_rank, color, CorePerNuma);
	if (MPI_Comm_split(_GLEXCOLL.Comm_intra, color, (_GLEXCOLL.global_rank), &(_GLEXCOLL.Comm_numa_intra)) != MPI_SUCCESS)
	{
		puts("Comm_intra MPI_Comm_split error");
		exit(0);
	}
	MPI_Comm_size(_GLEXCOLL.Comm_numa_intra, &(_GLEXCOLL.numa_intra_procn));
	MPI_Comm_rank(_GLEXCOLL.Comm_numa_intra, &(_GLEXCOLL.numa_intra_rank));

	// MPI_Barrier(MPI_COMM_WORLD);
	// if(_GLEXCOLL.global_rank == 0)
	// {
	// 	puts("==================================237=================================");
	// 	fflush(stdout);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	//根据每个进程分配numa inter communicator
	color = _GLEXCOLL.numa_intra_rank;
	if (MPI_Comm_split(_GLEXCOLL.Comm_Global, color, (_GLEXCOLL.global_rank), &(_GLEXCOLL.Comm_numa_inter)) != MPI_SUCCESS)
	{
		puts("Comm_numa_inter MPI_Comm_split error");
		exit(0);
	}
	MPI_Comm_size(_GLEXCOLL.Comm_numa_inter, &(_GLEXCOLL.numa_inter_procn));
	MPI_Comm_rank(_GLEXCOLL.Comm_numa_inter, &(_GLEXCOLL.numa_inter_rank));
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(_GLEXCOLL.global_rank == 0)
	// {
	// 	puts("=======================220=======================");
	// 	fflush(stdout);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// printf("host_name =%s global rank = %d,  numa_intra_rank = %d numa_inter_rank= %d\n",host_name,(_GLEXCOLL.global_rank),_GLEXCOLL.numa_intra_rank,_GLEXCOLL.numa_inter_rank);
	//  puts("check point");
	// numa 感知的多leader聚合通信卸载
	Comm_intra = (_GLEXCOLL.Comm_numa_intra);
	intra_procn = (_GLEXCOLL.numa_intra_procn);
	intra_rank = (_GLEXCOLL.numa_intra_rank);
	Comm_inter = (_GLEXCOLL.Comm_numa_inter);
	inter_rank = (_GLEXCOLL.numa_inter_rank);
	inter_procn = (_GLEXCOLL.numa_inter_procn);
	global_rank = (_GLEXCOLL.global_rank);
	global_procn = (_GLEXCOLL.global_procn);
	// printf("intra_procn =%d,inter_procn=%d\n ", intra_procn, inter_procn);

	//开始初始化线程有关的信息。
	// puts("check point225");
	pthread_mutex_init(&(_GLEXCOLL.mtx), NULL);
	sem_init(&(_GLEXCOLL.FinishREQ), 0, 1);
	sem_init(&(_GLEXCOLL.PostREQ), 0, 0);
	pthread_cond_init(&(_GLEXCOLL.cond), NULL);
	{
		//开始创建GLEX通信线程。
		_GLEXCOLL.th_flag = RUNNING;
		// puts("check point");
		pthread_create(&(_GLEXCOLL.comm_thread), NULL, comm_main_loop, NULL);
	}
	// puts("check point1");
	GLEX_ep_init();
	// puts("xe");
	for (int i = 0; i < GLEX_MELLOC_MAX; i++)
	{
		send_mhs[i] = (glex_mem_handle_t *)malloc(sizeof(glex_mem_handle_t) * inter_procn);
		recv_mhs[i] = (glex_mem_handle_t *)malloc(sizeof(glex_mem_handle_t) * inter_procn);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(_GLEXCOLL.global_rank == 0)
	{
		puts("=======================GLEXCOLL_Init初始化完成=======================");
		fflush(stdout);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// puts("end");
	//初始化allreduce
	// GLEXCOLL_InitAllreduce();
	//现在为每个节点上通信进程创建端点。
	//暂时为每个节点上主进程分配PIO端点。
}

void GLEXCOLL_Finalize()
{
	_GLEXCOLL.th_flag = EXIT;
	sem_post(&(_GLEXCOLL.PostREQ));
	pthread_mutex_destroy(&(_GLEXCOLL.mtx));
	pthread_cond_destroy(&(_GLEXCOLL.cond));
	//销毁通信线程
	pthread_join(_GLEXCOLL.comm_thread, NULL);
	GLEX_ep_destroy();
	for (int i = 0; i < GLEX_MELLOC_MAX; i++)
	{
		free(send_mhs[i]);
		free(recv_mhs[i]);
	}
}
