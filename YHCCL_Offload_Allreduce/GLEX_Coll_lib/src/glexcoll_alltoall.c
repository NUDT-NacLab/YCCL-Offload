#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "glexcoll.h"
#include "glexalltoall.h"
#include <sys/types.h>

#define PJT_ALLTOALL_DEBUG
//æ•´åž‹æ•°æ®åŽ‹ç¼©åº
#include "vint.h"
#include "vsimple.h"
#include "vp4.h"
#include "bitpack.h"
#include "eliasfano.h"

static char *bufS, *bufR;
static char *block_tmp;
//RDMAéœ€è¦ä½¿ç”¨çš„ç›¸å…³æ•°æ®
static struct glex_rdma_req rdma_req;
static struct glex_rdma_req *bad_rdma_req;
static glex_event_t *event;
static glex_mem_handle_t send_mh;
extern int leaderN;

extern void intra_memory_barrier_alltoall();
static int alltoall_global_rank;
static int alltoall_global_procn;
static int alltoall_inter_rank;
static int alltoall_inter_procn;
static int alltoall_intra_rank;
static MPI_Comm alltoall_comm_inter;
static MPI_Comm alltoall_comm_intra;
static int mesh_3d_count;
static int **mesh3dId_intra_mesh3dId_TO_a2a_grank;
static int host_id;
static int *alltoall_host_ids;
static int intra_Mesh3d_physical_id;
static int *intra_Mesh3d_physical_ids;
static int inter_Mesh3d_physical_id;
static int *inter_Mesh3d_physical_ids;
static int intra_rack_id;
static int inter_rack_id;
//TPDs17ç›¸å…³æ•°æ®
//æš‚æ—¶åªè€ƒè™‘æ¯èŠ‚ç‚¹å†…è¿›ç¨‹æ•°é‡ä¸^k
extern unsigned int get_coordinate_x_tpds(unsigned int in);
extern unsigned int get_coordinate_y_tpds(unsigned int in);
extern int Coordinate_x[64];
extern int Coordinate_y[64];
extern int Coordinate_numa_x[4];
extern int Coordinate_numa_y[4];

//å¤©æ²³3ç«‹æ–¹ä½“æ‹“æ‰‘æ„ŸçŸ¥è¡¨{cn0-cn7}->{<z,x,y>}
int physical_to_logical[8] = {0, 1, 3, 2, 7, 6, 4, 5};
int logical_to_physical[8] = {0, 1, 3, 2, 6, 7, 5, 4};
inline int mesh3Dlogical_flipZ(int in)
{
	return in ^ 4;
}
inline int mesh3Dlogical_flipX(int in)
{
	return in ^ 2;
}
inline int mesh3Dlogical_flipY(int in)
{
	return in ^ 1;
}
inline int mesh3Dlogical_flipZX(int in)
{
	return in ^ 6;
}
inline int mesh3Dlogical_flipZY(int in)
{
	return in ^ 5;
}
inline int mesh3Dlogical_flipXY(int in)
{
	return in ^ 3;
}

inline int mesh3Dlogical_flipXYZ(int in)
{
	return in ^ 7;
}

int mesh3D_physical_source_vec[8];
int mesh3D_physical_target_vec[8];
int Init_mesh3Dphysical_source_target(int physical)
{
	int logical_id = physical_to_logical[physical];
	mesh3D_physical_source_vec[0] = physical;
	mesh3D_physical_target_vec[0] = physical;

	mesh3D_physical_source_vec[1] = logical_to_physical[mesh3Dlogical_flipZ(logical_id)];
	mesh3D_physical_target_vec[1] = logical_to_physical[mesh3Dlogical_flipZ(logical_id)];

	mesh3D_physical_source_vec[2] = logical_to_physical[mesh3Dlogical_flipX(logical_id)];
	mesh3D_physical_target_vec[2] = logical_to_physical[mesh3Dlogical_flipX(logical_id)];

	mesh3D_physical_source_vec[3] = logical_to_physical[mesh3Dlogical_flipY(logical_id)];
	mesh3D_physical_target_vec[3] = logical_to_physical[mesh3Dlogical_flipY(logical_id)];

	mesh3D_physical_source_vec[4] = logical_to_physical[mesh3Dlogical_flipXYZ(logical_id)];
	mesh3D_physical_target_vec[4] = logical_to_physical[mesh3Dlogical_flipXYZ(logical_id)];

	mesh3D_physical_source_vec[5] = logical_to_physical[mesh3Dlogical_flipZX(logical_id)];
	mesh3D_physical_target_vec[5] = logical_to_physical[mesh3Dlogical_flipZX(logical_id)];

	mesh3D_physical_source_vec[6] = logical_to_physical[mesh3Dlogical_flipZY(logical_id)];
	mesh3D_physical_target_vec[6] = logical_to_physical[mesh3Dlogical_flipZY(logical_id)];

	mesh3D_physical_source_vec[7] = logical_to_physical[mesh3Dlogical_flipXY(logical_id)];
	mesh3D_physical_target_vec[7] = logical_to_physical[mesh3Dlogical_flipXY(logical_id)];
	if (alltoall_global_rank == 0)
		// if (0)
		if (intra_Mesh3d_physical_id == 1)
		{
			printf("ç«‹æ–¹ä½“å†…id = %d ä¸»æœºid=%d ä¸»æœºid mod 8=%d\n", intra_Mesh3d_physical_id, host_id, host_id % 8);
			for (int i = 0; i < 8; i++)
			{
				printf("source:cn%d target:cn%d\n", mesh3D_physical_source_vec[i], mesh3D_physical_target_vec[i]);
			}
		}
}

void Mesh3D_aware()
{
	//æ­¤å‡½æ•°ç”±æ¯ä¸ªèŠ‚ç‚¹çš„leadersè°ƒç”¨å³å¯
	//é€šä¿¡åŸŸä¸ºalltoall_comm_inter;
	int host_id;
	sscanf(host_name, "cn%d", &host_id);
	intra_Mesh3d_physical_id = host_id % 8;
	inter_Mesh3d_physical_id = host_id / 8;

	intra_rack_id = host_id % 256;
	inter_rack_id = host_id / 256;

	//ç¬¬ä¸€æ­¥ç”±leaderèŠ‚ç‚¹æ”¶é›†æ¯ä¸ªèŠ‚ç‚¹æ‰€åœ¨çš„3d meshåºå·
	inter_Mesh3d_physical_ids = (int *)malloc(sizeof(int) * alltoall_inter_procn);
	MPI_Allgather(&inter_Mesh3d_physical_id, 1, MPI_INT,
				  inter_Mesh3d_physical_ids, 1, MPI_INT, alltoall_comm_inter);

	intra_Mesh3d_physical_ids = (int *)malloc(sizeof(int) * alltoall_inter_procn);
	MPI_Allgather(&intra_Mesh3d_physical_id, 1, MPI_INT,
				  intra_Mesh3d_physical_ids, 1, MPI_INT, alltoall_comm_inter);

	//ç¬¬äºŒæ­¥ç»Ÿè®¡æ•´ä¸ªé€šä¿¡åŸŸä¸­Mesh3dçš„æ•°é‡
	mesh_3d_count = 1;
	for (int i = 1; i < alltoall_inter_procn; i++)
	{
		if (inter_Mesh3d_physical_ids[i - 1] != inter_Mesh3d_physical_ids[i])
			mesh_3d_count++;
	}
	//ç¬¬ä¸‰æ­¥åˆ†é…ç«‹æ–¹ä½“æ•°é‡ä¸ªäºŒç»´æ•°ç»„ã€
	mesh3dId_intra_mesh3dId_TO_a2a_grank = (int *)malloc(sizeof(int *) * mesh_3d_count);
	for (int i = 0; i < mesh_3d_count; i++)
	{
		mesh3dId_intra_mesh3dId_TO_a2a_grank[i] = (int *)malloc(sizeof(int) * 9);
		for (int j = 0; j < mesh3dId_intra_mesh3dId_TO_a2a_grank; j++)
		{
			mesh3dId_intra_mesh3dId_TO_a2a_grank[i][j] = -1;
		}
	}
	//ç¬¬å››æ­å»ºç«‹äºŒç»´æ˜ å°„è¡¨ï¼Œ[ç«‹æ–¹ä½“id][intra ç‰©ç†id]-> [a2a inter ranks]
	// {
	//     volatile int *pjt = (int *)malloc(sizeof(int));// = 0;
	//     *pjt = 0;
	//     pid_t pid = getpid();
	//     printf("%s %d\n", host_name, pid);
	//     while (*pjt == 0);
	// }
	//æ•°ç»„[èŠ¯ç‰‡é€»è¾‘id][ç«‹æ–¹ä½“é€»è¾‘id][ç«‹æ–¹ä½“å†…cn%8id 0-7]=mpirank
	{
		int j = 0, count = 1;

		int intra_phyid = intra_Mesh3d_physical_ids[0];
		mesh3dId_intra_mesh3dId_TO_a2a_grank[0][intra_phyid] = 0;
		for (int i = 1; i < alltoall_inter_procn; i++)
		{
			intra_phyid = intra_Mesh3d_physical_ids[0];
			if (inter_Mesh3d_physical_ids[i - 1] != inter_Mesh3d_physical_ids[i])
			{
				j++;
				mesh3dId_intra_mesh3dId_TO_a2a_grank[j][intra_phyid] = i;
			}
			else
			{
				mesh3dId_intra_mesh3dId_TO_a2a_grank[j][intra_phyid] = i;
			}
		}

		{
			//éªŒè¯æ˜ å°„è¡
			int start = 0;
			for (int row = 0; row < mesh_3d_count; row++)
			{
				int coln = mesh3dId_intra_mesh3dId_TO_a2a_grank[row][0];
				for (int i = 1; i <= coln; i++)
				{
					if (mesh3dId_intra_mesh3dId_TO_a2a_grank[row][i] != start)
					{
						printf("ç«‹æ–¹ä½“æ˜ å°„è¡¨æž„é€ é”™è¯meshindex=%d i=%d %d != %d\n", row, i, mesh3dId_intra_mesh3dId_TO_a2a_grank[row][i], start);
						exit(0);
					}
					start++;
				}
			}
		}
	}
	//ç¬¬å…­æ­å»ºç«‹ç«‹æ–¹ä½“å†…id intra_Mesh3d_physical_id åˆalltoall_inter_rankçš„æ˜ å°„è¡¨
	//ä¸»è¦ç”¨äºŽ ç«‹æ–¹ä½“å†…alltoallè½¬ç½®
	Init_mesh3Dphysical_source_target(intra_Mesh3d_physical_id);
}

int ***TopologyMap;
int chip_id;
int *chip_ids;
int mesh3d_id;
int *alltoall_send_order;
int *alltoall_recv_order;
int *alltoall_stepVec_send;
int *alltoall_stepVec_recv;
int alltoall_step_count;
void Topology_aware()
{
	// void Mesh3D_aware();
	//ç½‘ç»œè¢«åˆ†ä¸ºç«‹æ–¹ä½“ï¼ŒèŠ¯ç‰‡ä¸¤å±
	//å­˜å‚¨åœ¨TopologyMap[èŠ¯ç‰‡ç´¢å¼•][ç«‹æ–¹ä½“ç´¢å¼•][0-8] = mpirank
	int PJT_DEBUG = 0;
	int host_id;
	sscanf(host_name, "cn%d", &host_id);
	//TH3 æ¯28ä¸ªèŠ‚ç‚¹å…±äº«ä¸€ä¸ªè·¯ç”±èŠ¯ç‰
	// chip_id = host_id / 128;
	//æ¯ä¸ªèŠ‚ç‚¹å…±ä¸€ä¸d mesh
	mesh3d_id = (host_id % 128) / 8;
	//meshå†…id
	int node_id = (host_id % 128) % 8;
	// if (alltoall_intra_rank == 0)
	// {
	//     printf("chip_id=%d mesh3d_id=%d node_id=%d\n", host_id/256, mesh3d_id, node_id);
	// }

	//ç¬¬ä¸€æ­¥æ”¶é›†æ‰€æœ‰host_id
	int *host_ids = (int *)malloc(sizeof(int) * alltoall_inter_procn);
	MPI_Allgather(&host_id, 1, MPI_INT,
				  host_ids, 1, MPI_INT, alltoall_comm_inter);
	chip_ids = (int *)malloc(sizeof(int) * alltoall_inter_procn);
	MPI_Barrier(alltoall_comm_inter);
	// if (alltoall_global_rank == 0)
	//     puts("step 2");
	//ç¬¬äºŒæ­¥ç»Ÿè®¡èŠ¯ç‰‡æ€»æ•°é‡
	int chip_count = 1;
	chip_ids[0] = 0;
	for (int i = 1; i < alltoall_inter_procn; i++)
	{
		int prev_chipid = host_ids[i - 1] / 128;
		int his_chipid = host_ids[i] / 128;
		if (prev_chipid != his_chipid)
			chip_count++;
		chip_ids[i] = chip_count - 1;
	}
	TopologyMap = (int ***)malloc(sizeof(int **) * chip_count);
	for (int i = 0; i < chip_count; i++)
	{
		//th3æ¯ä¸ªèŠ¯ç‰‡å†…æœ‰16ä¸ªç«‹æ–¹ä½“
		TopologyMap[i] = (int **)malloc(sizeof(int *) * 16);
		for (int j = 0; j < 16; j++)
		{
			TopologyMap[i][j] = (int *)malloc(sizeof(int) * 8);
			for (int k = 0; k < 8; k++)
			{
				TopologyMap[i][j][k] = -1;
			}
		}
	}
	MPI_Barrier(alltoall_comm_inter);
	// if (alltoall_global_rank == 0)
	//     puts("step 3");
	//ç¬¬ä¸‰æ­¥å°†æ¯ä¸ªèŠ‚ç‚¹çš„æ‹“æ‰‘ä¿¡æ¯è£…å…¥TopologyMapä¸
	for (int i = 0; i < alltoall_inter_procn; i++)
	{
		int his_chip_id = chip_ids[i];
		int his_mesh_id = (host_ids[i] % 128) / 8;
		int his_node_id = (host_ids[i] % 8);
		// printf("%d %d %d\n",his_chip_id,his_mesh_id,his_node_id);
		TopologyMap[his_chip_id][his_mesh_id][his_node_id] = i;
		// if (alltoall_intra_rank == 0)
		// {
		//     printf("alltoall_inter_rank=%d hostname=%d chipid=%d mesh_id=%d node_id=%d\n",
		//             alltoall_inter_rank,
		//             host_ids[i],
		//             his_chip_id,
		//             his_mesh_id,
		//             his_node_id);
		// }
	}
	MPI_Barrier(alltoall_comm_inter);
	// if (alltoall_global_rank == 0)
	//     puts("step 4");
	int index = 0;
	int my_chip_id = chip_ids[alltoall_inter_rank];
	int my_mesh3d_id = (host_id % 128) / 8;
	int my_node_id = host_id % 8;
	Init_mesh3Dphysical_source_target(my_node_id);

	alltoall_stepVec_send = (int *)malloc(chip_count * 16 * sizeof(int));
	alltoall_stepVec_recv = (int *)malloc(chip_count * 16 * sizeof(int));

	// if(0)
	{
		for (int chipshift = 0; chipshift < chip_count; chipshift++)
		{
			for (int meshshift = 0; meshshift < 16; meshshift++)
			{
				int target_chip_id = (my_chip_id + chipshift) % chip_count;
				int target_mesh_id = (my_mesh3d_id + meshshift) % 16;
				int *target_vec = TopologyMap[target_chip_id][target_mesh_id];
				alltoall_stepVec_send[alltoall_step_count] = 0;
				for (int i = 0; i < 8; i++)
				{
					int ii = mesh3D_physical_source_vec[i]; //(i + my_node_id) % 8;
					if (target_vec[ii] != -1)
					{
						alltoall_send_order[index] = target_vec[ii];
						// if (alltoall_intra_rank == 0 && alltoall_inter_rank == 7)
						// {
						//     printf("%d <=> %d\n", alltoall_inter_rank, alltoall_send_order[index]);
						// }
						index++;
						alltoall_stepVec_send[alltoall_step_count] += 1;
					}
				}
				alltoall_step_count++;
			}
		}
		alltoall_step_count = 0;
		index = 0;
		for (int chipshift = 0; chipshift < chip_count; chipshift++)
		{
			for (int meshshift = 0; meshshift < 16; meshshift++)
			{
				int target_chip_id = (chip_count + my_chip_id - chipshift) % chip_count;
				int target_mesh_id = (16 + my_mesh3d_id - meshshift) % 16;
				int *target_vec = TopologyMap[target_chip_id][target_mesh_id];
				alltoall_stepVec_recv[alltoall_step_count] = 0;
				for (int i = 0; i < 8; i++)
				{
					int ii = mesh3D_physical_source_vec[i]; //(i + my_node_id) % 8;
					if (target_vec[ii] != -1)
					{
						alltoall_recv_order[index] = target_vec[ii];
						// if (alltoall_intra_rank == 0 && alltoall_inter_rank == 7)
						// {
						//     printf("%d <=> %d\n", alltoall_inter_rank, alltoall_recv_order[index]);
						// }
						index++;
						alltoall_stepVec_recv[alltoall_step_count] += 1;
					}
				}
				alltoall_step_count++;
			}
		}
	}
	if (PJT_DEBUG)
		printf("alltoall_step_countxx = %d\n", alltoall_step_count);

	if (index != alltoall_inter_procn)
	{
		puts("error: é‡æŽ’åŽæ¶ˆæ¯æ€»é‡ä¸Šä¸åŽ»ï¼Ÿ");
		exit(0);
	}
}
glex_mem_handle_t *alltoall_mem_handle_vec;
glex_ep_addr_t *alltoall_ep_addrs;
long long BUFR_size = (1L << 28);
long long BUFS_size = (1L << 28);
int BUFR_block_size = (1 << 28);
int BUFS_block_size = (1 << 28);
static glex_mem_handle_t *alltoall_recv_mhs;
static glex_mem_handle_t **alltoall_recv_mhs_vec;
static glex_mem_handle_t *alltoall_send_mhs;
int BUFR_block_count;
int BUFS_block_count;
static void glexcoll_Init_RDMA(MPI_Comm comm)
{
	BUFR_block_count = 1; //(BUFR_size / BUFR_block_size);
	BUFS_block_count = 1; //(BUFS_size / BUFS_block_size);

	// leaderN = mmax(mmin(ppn/4,4),1);
	// printf("leaderN = %d\n",leaderN);
#ifdef PJT_NEW_VERSION
	if (alltoall_inter_procn > 1 && am_i_leader())
#else
	if (alltoall_inter_procn > 1 && alltoall_intra_rank % 4 == 0 && (alltoall_intra_rank >> 2) < leaderN)
#endif

	// if (alltoall_inter_procn > 1 && alltoall_intra_rank % 4 == 0 && (alltoall_intra_rank >> 2) < leaderN)
	{
		// alltoall_mem_handle_vec = (glex_mem_handle_t *)malloc(sizoef(glex_mem_handle_t) * alltoall_inter_procn);

		alltoall_ep_addrs = (glex_ep_addr_t *)malloc(sizeof(glex_ep_addr_t) * alltoall_inter_procn);
		// MPI_Allgather((viod *)&(_GLEXCOLL.ep), );
		{
			glex_ep_addr_t my_ep_addr;
			glex_get_ep_addr(_GLEXCOLL.ep, &my_ep_addr);
			MPI_Allgather((void *)&my_ep_addr, sizeof(my_ep_addr), MPI_CHAR, alltoall_ep_addrs, sizeof(my_ep_addr), MPI_CHAR, alltoall_comm_inter);
		}
		bufR = (char *)malloc(BUFR_size);
		bufS = (char *)malloc(BUFS_size);
		alltoall_recv_mhs = (glex_mem_handle_t *)malloc(sizeof(*alltoall_recv_mhs) * BUFR_block_count);
		alltoall_recv_mhs_vec = (glex_mem_handle_t **)malloc(sizeof(*alltoall_recv_mhs_vec) * BUFR_block_count);
		alltoall_send_mhs = (glex_mem_handle_t *)malloc(sizeof(*alltoall_send_mhs) * BUFS_block_count);

		// MPI_Barrier(alltoall_comm_inter);
		// 	if(inter_rank == 0 && intra_rank == 0) puts("419");
		// printf("global_rank = %d, %p\n",global_rank,_GLEXCOLL.ep);
		//将多个接收BLOCK缓冲区分块注册。
		for (int i = 0; i < BUFR_block_count; i++)
		{
			void *startp = bufR + i * BUFR_block_size;
			int ret = glex_register_mem(_GLEXCOLL.ep, startp, BUFR_block_size,
										GLEX_MEM_WRITE,
										&(alltoall_recv_mhs[i]));
			if (ret != GLEX_SUCCESS)
			{
				printf("_register_mem(), return: %d\n", ret);
				perror("");
				exit(1);
			}
			alltoall_recv_mhs_vec[i] = (glex_mem_handle_t *)malloc(sizeof(_GLEXCOLL.local_mh) * alltoall_inter_procn);
			MPI_Allgather(&(alltoall_recv_mhs[i]), sizeof(_GLEXCOLL.local_mh), MPI_CHAR,
						  alltoall_recv_mhs_vec[i], sizeof(_GLEXCOLL.local_mh), MPI_CHAR, alltoall_comm_inter);
			// MPI_Barrier(alltoall_comm_inter);
			// if(inter_rank == 0 && intra_rank == 0)puts("438");
		}
		_GLEXCOLL.local_mh = alltoall_recv_mhs[0];
		_GLEXCOLL.rmt_mhs = (glex_mem_handle_t *)malloc(sizeof(_GLEXCOLL.local_mh) * alltoall_inter_procn);
		MPI_Allgather((void *)&(_GLEXCOLL.local_mh), sizeof(_GLEXCOLL.local_mh), MPI_CHAR,
					  _GLEXCOLL.rmt_mhs, sizeof(_GLEXCOLL.local_mh), MPI_CHAR, alltoall_comm_inter);
		block_tmp = (char *)malloc(ppn * (1 << 20));

		// MPI_Barrier(alltoall_comm_inter);
		// 	if(inter_rank == 0 && intra_rank == 0)puts("442");
		//将多个发送BLOCK缓冲区注册
		for (int i = 0; i < BUFS_block_count; i++)
		{
			void *startp = bufS + i * BUFS_block_size;
			int ret = glex_register_mem(_GLEXCOLL.ep, startp, BUFS_block_size,
										GLEX_MEM_READ,
										&(alltoall_send_mhs[i]));
			if (ret != GLEX_SUCCESS)
			{
				printf("_register_mem(), return: %d\n", ret);
				perror("");
				exit(1);
			}
			// MPI_Barrier(alltoall_comm_inter);
			// if(inter_rank == 0 && intra_rank == 0) puts("461");
		}
		send_mh = alltoall_send_mhs[0];
		// // for(int i = 0;i<BUFR_block_count;)
		// int ret = glex_register_mem(_GLEXCOLL.ep, bufS, BUFSR_size,
		// 		GLEX_MEM_READ | GLEX_MEM_WRITE,
		// 		&send_mh);
		// if (ret != GLEX_SUCCESS)
		// {
		// 	printf("_register_mem(), return: %d\n", ret);
		// 	perror("");
		// 	exit(1);
		// }

		// printf("%d check Init RDMA]\n",alltoall_global_rank);
	}
	else
	{
		bufR = (char *)malloc(128);
		bufS = (char *)malloc(128);
	}
}

void glexcoll_register_alltoall_buffer_new(void *sendbuf, void *recvbuf, int size, MPI_Comm comm, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	if (ppn != 1)
	{
		puts("每节点多进程请用  glexcoll_register_alltoall_buffer");
		exit(0);
	}
	int procn;
	MPI_Comm_size(comm, &procn);
	if (size > 1 << 27)
	{
		puts("单缓冲区大小超过128MB，可能导致无法注册内存");
	}
	bufmh->alltoall_ep_addrs = (glex_ep_addr_t *)malloc(sizeof(glex_ep_addr_t) * procn);
	{
		glex_ep_addr_t my_ep_addr;
		glex_get_ep_addr(_GLEXCOLL.ep, &my_ep_addr);
		MPI_Allgather((void *)&my_ep_addr, sizeof(my_ep_addr), MPI_CHAR, bufmh->alltoall_ep_addrs, sizeof(my_ep_addr), MPI_CHAR, comm);
	}
	glex_mem_handle_t tmp;
	int ret = glex_register_mem(_GLEXCOLL.ep, recvbuf, size,
								GLEX_MEM_READ | GLEX_MEM_WRITE,
								&(tmp));
	ret = glex_register_mem(_GLEXCOLL.ep, sendbuf, size,
							GLEX_MEM_READ | GLEX_MEM_WRITE,
							&(bufmh->sendmh));
	if (ret != GLEX_SUCCESS)
	{
		printf("_register_mem(), return: %d\n", ret);
		perror("");
		exit(1);
	}
	bufmh->alltoall_mem_handle_vec = (glex_mem_handle_t *)malloc(sizeof(tmp) * procn);
	MPI_Allgather((void *)&(tmp), sizeof(tmp), MPI_CHAR,
				  bufmh->alltoall_mem_handle_vec, sizeof(tmp), MPI_CHAR, comm);
	bufmh->sendbuf = sendbuf;
	bufmh->recvbuf = recvbuf;
}

static void glexcoll_Destroy_RDMA()
{
	// if (alltoall_inter_procn > 1 && alltoall_intra_rank % 4 == 0 && (alltoall_intra_rank >> 2) < leaderN)
#ifdef PJT_NEW_VERSION
	if (alltoall_inter_procn > 1 && am_i_leader())
#else
	if (alltoall_inter_procn > 1 && alltoall_intra_rank % 4 == 0 && (alltoall_intra_rank >> 2) < leaderN)
#endif
	{
		free(bufR);
		free(bufS);
	}
}

extern void glexcoll_init_alltoall_shared_memory_buffer();
extern void TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(int shift, int elem_size);
extern void TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(int shift, int elem_size, struct GLEXCOLL_a2a_bufmh *bufmh);
//åˆå§‹åŒ–è¿‡ç¨‹ï¼Œæ³¨æ„commä¸­å¿…é¡»å®Œæ•´å¡«æ»¡æ¯ä¸€ä¸ªèŠ‚ç‚
void glexcoll_InitAlltoall_new(MPI_Comm comm)
{
	MPI_Comm_rank(comm, &alltoall_global_rank);
	MPI_Comm_size(comm, &alltoall_global_procn);
	alltoall_intra_rank = intra_rank;

	int color = alltoall_intra_rank;
	MPI_Comm_split(comm, color, alltoall_global_rank, &(alltoall_comm_inter));
	MPI_Comm_rank(alltoall_comm_inter, &alltoall_inter_rank);
	MPI_Comm_size(alltoall_comm_inter, &alltoall_inter_procn);
	// printf("alltoall_global_rank=%d alltoall_global_procn=%d alltoall_inter_rank=%d alltoall_inter_procn=%d alltoall_intra_rank=%d\n",
	//         alltoall_global_rank, alltoall_global_procn, alltoall_inter_rank, alltoall_inter_procn, alltoall_intra_rank);
	// glexcoll_init_alltoall_shared_memory_buffer();
	if (ppn > 1)
		glexcoll_Init_RDMA(comm);

	//åˆå§‹åŒ–NUMA-Coæ‰€éœ€è¦çš„è½¬ç½®é¡ºåº
	{
		int CorePnuma = 4;
		int start_index = (intra_rank % CorePnuma) * 4;
		int end_index = start_index + 3;
		for (int i = start_index; i <= end_index; i++)
		{
			//åˆå§‹åŒ–è½¬ç½®ç›®æ ‡çš„åæ ‡
			Coordinate_numa_x[i - start_index] = get_coordinate_x_tpds(i);
			Coordinate_numa_y[i - start_index] = get_coordinate_y_tpds(i);
			// if(intra_rank == 7)
			// {
			//     printf("global_rank=%d, index = %d %d %d \n",global_rank,i,Coordinate_numa_x[i-start_index],Coordinate_numa_y   [i-start_index]);
			// }
		}
		// exit(0);
	}
	//åˆå§‹åŒ–ç«‹æ–¹ä½“æ„ŸçŸ¥alltoallæ‰€éœ€çš„æ•°æ®ç»“æž
	{
		alltoall_send_order = (int *)malloc(sizeof(int) * alltoall_inter_procn);
		alltoall_recv_order = (int *)malloc(sizeof(int) * alltoall_inter_procn);
		// printf("global_rank = %d",global_rank);
		if (alltoall_intra_rank == 0)
		{
			// Mesh3D_aware();
			Topology_aware();
		}
		MPI_Bcast(alltoall_send_order, alltoall_inter_procn, MPI_INT, 0, Comm_intra);
		MPI_Bcast(alltoall_recv_order, alltoall_inter_procn, MPI_INT, 0, Comm_intra);
	}
}

void glexcoll_intra_rma_gather(void *bufrecv, int send_shift, int msgsz)
{
	char *p = bufrecv;
	for (int s = 0; s < intra_procn; s++)
	{
		volatile char *source = get_is_senddata_buffer(s) + send_shift;
		for (int i = 0; i < msgsz; i++)
		{
			p[i] = source[i];
		}
		p += msgsz;
	}
}
void glexcoll_intra_rma_scatter(void *bufsend, int recv_shift, int msgsz)
{
	// {
	//     volatile int pjt = 0;
	//     pid_t pid = getpid();
	//     printf("%s %d\n", host_name, pid);
	//     while (pjt == 0) ;
	// }

	char *p = bufsend;
	for (int s = 0; s < intra_procn; s++)
	{
		volatile char *target = get_is_recvdata_buffer(s) + recv_shift;
		for (int i = 0; i < msgsz; i++)
		{
			target[i] = p[i];
		}
		p += msgsz;
	}
}

void MATRIX_transform_local(char *buf, int n, int block_size)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			char *p = buf + (i * n + j) * block_size;
			char *q = buf + (j * n + i) * block_size;
			memcpy(block_tmp, p, block_size);
			memcpy(p, q, block_size);
			memcpy(q, block_tmp, block_size);
			// printf("global_rank=%d *q = %d\n",global_rank, *(int*)q);
		}
	}
}
void glexcoll_multileader_pipeline_rdma_alltoall(void *sendbuf,
												 int sendsize,
												 void *recvbuf,
												 int recvsize,
												 MPI_Comm comm)
{

	extern int Intra_Gather_Scatter_ON;
	// puts("129");
	extern int num_of_ongoing_msg;
	// MPI_Alltoall(sendbuf, sendsize, MPI_CHAR, recvbuf, recvsize, MPI_CHAR, comm);
	// MPI_Barrier(comm);
	//ç¬¬ä¸€æ­¥å…ˆè¿›è¡ŒèŠ‚ç‚¹å†…alltoallè½¬ç½®
	int block_size = sendsize * ppn;
	int shift = ppn * sendsize * alltoall_inter_rank;
	// int *p = get_is_senddata_buffer(intra_rank);
	// {
	//     //æ‰“å°è¾“å‡ºèŠ‚ç‚¹å†…è½¬ç½®ç»“æžœã€
	//     //2èŠ‚ç‚¹xppnæ ¸å¿ƒã€
	//     //æ£€æŸ¥æ­£ç¡®ã€‚z
	//     if (global_rank >= 8)
	//     {
	//         puts("----------------------------");
	//         for (int i = 8; i < 9; i++)
	//         {
	//             printf("send : %d == %d\n", ((int *)(sendbuf))[i], p[i]);
	//         }
	//     }
	// }
	// TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(shift, sendsize);
	TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(shift, sendsize);
	// MPI_Barrier(comm);
	// if (alltoall_global_rank == 0)
	//     puts("check a start multi leader a2a");
	// MPI_Barrier(comm);
	// TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(shift, sendsize);
	// {
	//     //æ‰“å°è¾“å‡ºèŠ‚ç‚¹å†…è½¬ç½®ç»“æžœã€
	//     //2èŠ‚ç‚¹xppnæ ¸å¿ƒã€
	//     //æ£€æŸ¥æ­£ç¡®ã€‚z
	//     if (global_rank == 8)
	//     {
	//         puts("----------------------------");
	//         for (int i = 8; i < 16; i++)
	//         {
	//             printf("recv %d ", ((int *)(recvbuf))[i]);
	//         }
	//         puts("\n----------------------------");
	//     }
	// }

	// TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(shift, sendsize);
	//ç¬¬äºŒæ­¥è¿›è¡ŒèŠ‚ç‚¹é—´alltoallè½¬ç½®
	//é¦–å…ˆè¿›è¡ŒGatheræ“ä½œï¼Œå°†æ•°æ®ç”±å¤šä¸ªLeaderå„è‡ª Gatheråˆ°è‡ªå·±çš„BufS
	int BufScount = 0;
	{
		char *p = bufS;
		for (int shift = 1; shift < alltoall_inter_procn; ++shift)
		{
			int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
			int leader = (shift % leaderN) << 2;
			int source = sendbuf + block_id * block_size;
			if (alltoall_intra_rank == leader)
			{
				if (Intra_Gather_Scatter_ON)
					glexcoll_intra_rma_gather(p, block_id * block_size, block_size);
				{
					//æ‰“å°gatherè¾“å‡ºæ•°æ®
					// if (global_rank == 12)
					// {
					//     puts("\n-------------------------------");
					//     for (int i = 0; i < ppn * ppn; i++)
					//     {
					//         printf("gather result: global_rank=%d %d\n", global_rank, ((int *)p)[i]);
					//     }
					// }
				}
				//æ•°æ®gatherå®ŒæˆåŽç«‹é©¬å‘é€å‡ºåŽ
				{
					// {
					//     volatile int pjt = 0;
					//     printf("%s %d\n",host_name,alltoall_intra_rank);
					//     while(pjt == 0) ;
					// }
					int target = (alltoall_inter_rank + shift) % alltoall_inter_procn;
					rdma_req.rmt_ep_addr.v = alltoall_ep_addrs[target].v; //_GLEXCOLL.ep_addrs[target].v;
					rdma_req.local_mh.v = send_mh.v;
					rdma_req.local_offset = (p - bufS);
					rdma_req.len = block_size * ppn;
					rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
					rdma_req.rmt_offset = (BufScount)*ppn * block_size;
					rdma_req.type = GLEX_RDMA_TYPE_PUT;
					rdma_req.rmt_evt.cookie[0] = 2 + alltoall_inter_rank;
					rdma_req.rmt_evt.cookie[1] = 2 + BufScount;
					rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
					// if (BufScount % num_of_ongoing_msg == 1)
					//     rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
					// else
					rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
					rdma_req.next = NULL;
					int ret;
					while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
					{
					}
					if (ret != GLEX_SUCCESS)
					{
						if (ret == GLEX_INVALID_PARAM)
							printf("%d, _rdma() 非法参数", global_rank);
						printf("_rdma(), return: %d\n", ret);
						exit(1);
					}
					BufScount++;
				}
				//å‘é€æŒ‡é’ˆæ›´æ–
				p += block_size * ppn;
				// printf("%d check send\n",alltoall_global_rank);
			}
		}
		{

			//leaderèŠ‚ç‚¹å¿…é¡»ç­‰å¾…æ‰€æœ‰æ¶ˆæ¯åˆ°è¾
			for (int c = 0; c < BufScount; ++c)
			{
				while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
					;

				// {
				//     volatile int pjt = 0;
				//     printf("check recv %s %d\n",host_name,alltoall_intra_rank);
				//     while(pjt == 0) ;
				// }

				int source_inter_rank = event->cookie[0] - 1;
				int buf_recv_shift = event->cookie[1] - 1;
				// printf("%d recv a msg from source_inter_rank=%d\n",global_rank,source_inter_rank);
				//leaderèŠ‚ç‚¹è´Ÿè´£å°†æ¶ˆæ¯è½¬ç½
				int *source = bufR + buf_recv_shift * block_size * ppn;
				// if(global_rank == 12)
				// {
				//     for (int i = 0; i < ppn*ppn; i++)
				//     {
				//         printf("è½¬ç½®å‰ç»“æžœglobal_rank=%d %d\n", global_rank, ((int *)source)[i]);
				//     }
				//     puts("---------------------");
				// }
				// printf("start MATRIX_transform_local %d\n",global_rank);
				extern int MATRIX_Tanfform_ON;
				if (MATRIX_Tanfform_ON)
				{
					MATRIX_transform_local(source, ppn, sendsize);
					// if(global_rank == 4)
					// {
					//     for (int i = 0; i < ppn*ppn; i++)
					//     {
					//         printf("è½¬ç½®åŽç»“æžœglobal_rank=%d %d\n", global_rank, ((int *)source)[i]);
					//     }
					// }
				}
				// printf("start glexcoll_intra_rma_scatter %d BufScount=%d\n",global_rank,BufScount);
				//è½¬ç½®ä¹‹åŽï¼Œå°†æ¶ˆæ¯Scatteråˆ°å¯¹åº”å¤„
				// if(Intra_Gather_Scatter_ON)
				{
					glexcoll_intra_rma_scatter(source, source_inter_rank * block_size, block_size);
				}
			}
		}
		if (BufScount > 0)
			glex_discard_probed_event(_GLEXCOLL.ep);
		// MPI_Barrier(comm);
		// if(alltoall_global_rank == 0) puts("check all rdma");
	}

	// printf("start barrier %d\n", global_rank);
	intra_memory_barrier_alltoall();

	// for (int i = 0; i < alltoall_global_procn; i++)
	// {
	//         printf("%d recv%d\n", global_rank, ((int *)recvbuf)[i]);
	// }
	// intra_memory_barrier_alltoall();
	// exit(0);
	//æŒ‰ç…§Linear shift exchangeçš„æ–¹å¼è¿›è¡Œgather
	// MPI_Barrier(comm);
}

extern void *Get_sendbuf(struct GLEXCOLL_a2a_bufmh *mh, unsigned long long i);
extern void *Get_recvbuf(struct GLEXCOLL_a2a_bufmh *mh, unsigned long long i);

void glexcoll_intra_rma_gather_with_mh(void *bufrecv, int send_shift, int msgsz, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	char *p = bufrecv;
	for (int s = 0; s < intra_procn; s++)
	{
		volatile char *source = Get_sendbuf(bufmh, s) + send_shift;
		// printf("%d gather %d's data send_shift=%d\n", alltoall_global_rank,s,send_shift);

		for (int i = 0; i < msgsz; i++)
		{
			p[i] = source[i];
		}
		p += msgsz;
	}
}

void glexcoll_intra_rma_scatter_with_mh(void *bufsend, int recv_shift, int msgsz, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	char *p = bufsend;
	for (int s = 0; s < intra_procn; s++)
	{
		volatile char *target = Get_recvbuf(bufmh, s) + recv_shift;
		for (int i = 0; i < msgsz; i++)
		{
			target[i] = p[i];
		}
		p += msgsz;
	}
}
void glexcoll_intra_rma_gather_NUMA_Aware_with_mh(void *bufrecv, int send_shift, int msgsz, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	char *p = bufrecv;
	for (int s = 0; s < intra_procn; s++)
	{
		int source_rank = (alltoall_intra_rank + s) % intra_procn;
		volatile char *source = Get_sendbuf(bufmh, source_rank) + send_shift;
		// printf("%d gather %d's data send_shift=%d\n", alltoall_global_rank,s,send_shift);
		char *p_t = p + source_rank * msgsz;
		for (int i = 0; i < msgsz; i++)
		{
			p_t[i] = source[i];
		}
	}
}
void glexcoll_intra_rma_scatter_NUMA_Aware_with_mh(void *bufsend, int recv_shift, int msgsz, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	char *p = bufsend;
	for (int s = 0; s < intra_procn; s++)
	{
		int target_rank = (alltoall_intra_rank + s) % intra_procn;
		volatile char *target = Get_recvbuf(bufmh, target_rank) + recv_shift;
		char *p_t = p + target_rank * msgsz;
		for (int i = 0; i < msgsz; i++)
		{
			target[i] = p_t[i];
		}
	}
}

extern int mmin(int a, int b);
extern int mmax(int a, int b);

void glexcoll_ONMPML_a2a(int sendsize,
						 struct GLEXCOLL_a2a_bufmh *bufmh)
{
	//此处代码发生段错误
	void *sendbuf = Get_sendbuf(bufmh, intra_rank);
	void *recvbuf = Get_recvbuf(bufmh, intra_rank);

	extern int Intra_Gather_Scatter_ON;
	extern int num_of_ongoing_msg;
	//ç¬¬ä¸€æ­¥å…ˆè¿›è¡ŒèŠ‚ç‚¹å†…alltoallè½¬ç½®
	int block_size = sendsize * ppn;
	int shift = ppn * sendsize * alltoall_inter_rank;
	if (ppn != 1)
	{
		TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(shift, sendsize, bufmh);
	}
	else
	{
		memcpy(recvbuf + alltoall_inter_rank * sendsize, sendbuf + alltoall_inter_rank * sendsize, sendsize);
	}

	// printf("sendbuf[myself]=%d, recvbuf[myself=%d]=%d\n", ((int *)sendbuf)[alltoall_inter_rank], alltoall_inter_rank, ((int *)recvbuf)[alltoall_inter_rank]);
	int inter_node_size = block_size * ppn;
	int steps_per_rounds = leaderN * ((BUFS_block_size) / inter_node_size);

	for (unsigned int h = 1; h < alltoall_inter_procn; h += steps_per_rounds)
	{
		// MPI_Barrier(MPI_COMM_WORLD);    if(alltoall_global_rank == 0) printf("check 721\n");

		int start = h;
		int end = mmin(alltoall_inter_procn, start + steps_per_rounds);
		char *p = bufS;
		int BufScount = 0;
		// if(global_rank == 0)
		// 	printf("ONMPMP start= %d end = %d\n",start,end);
		for (int shift = start; shift < end; ++shift)
		{
			// MPI_Barrier(MPI_COMM_WORLD);
			int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
			// if(global_rank == 4 || global_rank == 25)
			// 	printf("rank = %d 923 shift mod leaderN = %d\n",global_rank,shift % leaderN);
			int at = shift % leaderN;
			// if(global_rank == 4 || global_rank == 25)
			// 	printf("rank = %d 926 shift mod leaderN = %d\n",global_rank,shift % leaderN);
			int leader = get_ith_leader(at);
			// if(global_rank == 4 || global_rank == 25)
			// 	printf("rank = %d 929 shift mod leaderN = %d\n",global_rank,shift % leaderN);
			if (alltoall_intra_rank == leader)
			{
				// puts("check 930");
				if (p - bufS >= BUFS_block_size)
				{
					printf("ONMPMP rdma发送缓冲区溢出,global rank = %d shift=%d leader =%d\n", global_rank, shift, leader);
					// exit(0);
				}

				glexcoll_intra_rma_gather_NUMA_Aware_with_mh(p, block_id * block_size, block_size, bufmh);
				// if(alltoall_global_rank == 30) puts("check 932");
				{
					int target = (alltoall_inter_rank + shift) % alltoall_inter_procn;
					// if(alltoall_global_rank == 30) puts("check 943");

					rdma_req.rmt_ep_addr.v = alltoall_ep_addrs[target].v; //_GLEXCOLL.ep_addrs[target].v;
					// printf("rank = %d check 760\n", alltoall_global_rank);
					rdma_req.local_mh.v = send_mh.v;
					rdma_req.local_offset = (p - bufS);
					rdma_req.len = block_size * ppn;
					rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
					rdma_req.rmt_offset = BufScount * ppn * block_size;
					rdma_req.type = GLEX_RDMA_TYPE_PUT;
					// if(alltoall_global_rank == 30) puts("check 953");
					rdma_req.rmt_evt.cookie[0] = 1 + alltoall_inter_rank;
					rdma_req.rmt_evt.cookie[1] = 1 + BufScount;
					rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
					if (BufScount % num_of_ongoing_msg == 1)
						rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
					else
						rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
					rdma_req.next = NULL;
					int ret;
					while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
					{
					}
					if (ret != GLEX_SUCCESS)
					{
						if (ret == GLEX_INVALID_PARAM)
							printf("%d, _rdma() 非法参数", global_rank);
						printf("_rdma(), return: %d\n", ret);
						exit(1);
					}
					BufScount++;
				}
				// if(alltoall_global_rank == 30) puts("check 974");

				p += block_size * ppn;
			}
		}

		// if(global_rank == 4)
		// 	printf("rank = %d check 758\n",alltoall_global_rank);
		for (int c = 0; c < BufScount; ++c)
		{
			while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
				;

			// if(alltoall_global_rank == 30) puts("check 988");
			int source_inter_rank = event->cookie[0] - 1;
			int buf_recv_shift = event->cookie[1] - 1;
			//leaderèŠ‚ç‚¹è´Ÿè´£å°†æ¶ˆæ¯è½¬ç½
			int *source = bufR + buf_recv_shift * block_size * ppn;
			extern int MATRIX_Tanfform_ON;
			if (MATRIX_Tanfform_ON)
			{
				MATRIX_transform_local(source, ppn, sendsize);
			}
			//è½¬ç½®ä¹‹åŽï¼Œå°†æ¶ˆæ¯Scatteråˆ°å¯¹åº”å¤„
			// if(global_rank == 25 || global_rank == 1)
			// printf("%d recv from %d buf_recv_shift=%d\n",alltoall_inter_rank,source_inter_rank,buf_recv_shift);
			{
				glexcoll_intra_rma_scatter_NUMA_Aware_with_mh(source, source_inter_rank * block_size, block_size, bufmh);
				// glexcoll_intra_rma_scatter_with_mh(source, source_inter_rank * block_size, block_size, bufmh);
			}
		}

		if (am_i_leader())
		{
			glex_discard_probed_event(_GLEXCOLL.ep);
			// if(inter_rank == 0)
			// 	puts("1368");
			if (h + steps_per_rounds < alltoall_inter_procn)
				MPI_Barrier(alltoall_comm_inter);
			// if(inter_rank == 0)
			// 	puts("1370");
		}
	}

	// if(global_rank == 25 || global_rank == 4)
	// 	printf("rank = %d check 758\n",alltoall_global_rank);
	intra_memory_barrier_alltoall();

	// MPI_Barrier(MPI_COMM_WORLD);    if(alltoall_global_rank == 4) puts("check 1006");
	// intra_memory_barrier_alltoall();
	// printf("sendbuf[myself]=%d, recvbuf[myself=%d]=%d\n", ((int *)sendbuf)[alltoall_inter_rank], alltoall_inter_rank, ((int *)recvbuf)[alltoall_inter_rank]);
}

void glexcoll_SONMPML_a2a_with_mh(int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	void *sendbuf = Get_sendbuf(bufmh, intra_rank);
	void *recvbuf = Get_recvbuf(bufmh, intra_rank);

	// puts("check 710");
	// MPI_Barrier(MPI_COMM_WORLD);
	// puts("129");
	extern int Intra_Gather_Scatter_ON;
	extern int num_of_ongoing_msg;
	//ç¬¬ä¸€æ­¥å…ˆè¿›è¡ŒèŠ‚ç‚¹å†…alltoallè½¬ç½®
	int block_size = sendsize * ppn;
	int shift = ppn * sendsize * alltoall_inter_rank;
	if (ppn != 1)
	{
		TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(shift, sendsize, bufmh);
	}
	else
	{
		memcpy(recvbuf + alltoall_inter_rank * sendsize, sendbuf + alltoall_inter_rank * sendsize, sendsize);
	}

	int slice_size = 8192;
	int slice_count = block_size * ppn / slice_size;

	int BufScount = 0;
	{
		char *p = bufS;
		for (int shift = 1; shift < alltoall_inter_procn; ++shift)
		{
			int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
			int leader = get_ith_leader(shift % leaderN);
			if (alltoall_intra_rank == leader)
			{
				glexcoll_intra_rma_gather_NUMA_Aware_with_mh(p, block_id * block_size, block_size, bufmh);
				{
					//æ‰“å°gatherè¾“å‡ºæ•°æ®
					// if (alltoall_inter_rank == 0)
					// {
					//     puts("\n-------------------------------");
					//     for (int i = 0; i < ppn * ppn; i++)
					//     {
					//         printf("gather result: global_rank=%d %d\n", global_rank, ((int *)p)[i]);
					//     }
					// }
				}
				// MPI_Barrier(alltoall_comm_inter);
				//æ•°æ®gatherå®ŒæˆåŽç«‹é©¬å‘é€å‡ºåŽ
				int start_shift = 0;
				for (int h = 0; h < slice_count; h++)
				{
					int len = 0;
					if (h == slice_count - 1)
						len = block_size * ppn - start_shift;
					else
					{
						len = slice_size;
					}
					int target = (alltoall_inter_rank + shift) % alltoall_inter_procn;
					// printf("rank = %d check 758\n",alltoall_global_rank);
					rdma_req.rmt_ep_addr.v = alltoall_ep_addrs[target].v; //_GLEXCOLL.ep_addrs[target].v;
					// printf("rank = %d check 760\n", alltoall_global_rank);
					rdma_req.local_mh.v = send_mh.v;
					rdma_req.local_offset = start_shift + (p - bufS);
					rdma_req.len = len;
					rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
					rdma_req.rmt_offset = start_shift + BufScount * ppn * block_size;
					rdma_req.type = GLEX_RDMA_TYPE_PUT;
					rdma_req.rmt_evt.cookie[0] = 1 + alltoall_inter_rank;
					rdma_req.rmt_evt.cookie[1] = 1 + BufScount;
					rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
					rdma_req.flag = NULL;
					if (h == slice_count - 1)
					{
						rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
					}
					rdma_req.next = NULL;
					int ret;
					while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
					{
					}
					if (ret != GLEX_SUCCESS)
					{
						if (ret == GLEX_INVALID_PARAM)
							printf("%d, _rdma() 非法参数", global_rank);
						printf("_rdma(), return: %d\n", ret);
						exit(1);
					}
					start_shift += h * slice_size;
				}
				BufScount++;
				p += block_size * ppn;
			}
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);    puts("check 788");

	{
		//leaderèŠ‚ç‚¹å¿…é¡»ç­‰å¾…æ‰€æœ‰æ¶ˆæ¯åˆ°è¾
		for (int c = 0; c < BufScount; ++c)
		{
			while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
				;

			int source_inter_rank = event->cookie[0] - 1;
			int buf_recv_shift = event->cookie[1] - 1;
			//leaderèŠ‚ç‚¹è´Ÿè´£å°†æ¶ˆæ¯è½¬ç½
			int *source = bufR + buf_recv_shift * block_size * ppn;
			extern int MATRIX_Tanfform_ON;
			if (MATRIX_Tanfform_ON)
			{
				MATRIX_transform_local(source, ppn, sendsize);
			}
			//è½¬ç½®ä¹‹åŽï¼Œå°†æ¶ˆæ¯Scatteråˆ°å¯¹åº”å¤„
			// printf("%d recv from %d buf_recv_shift=%d\n",alltoall_inter_rank,source_inter_rank,buf_recv_shift);
			{
				glexcoll_intra_rma_scatter_NUMA_Aware_with_mh(source, source_inter_rank * block_size, block_size, bufmh);
				// glexcoll_intra_rma_scatter_with_mh(source, source_inter_rank * block_size, block_size, bufmh);
			}
		}
	}
	if (BufScount > 0)
		glex_discard_probed_event(_GLEXCOLL.ep);

	intra_memory_barrier_alltoall();
	if (alltoall_intra_rank == 0)
		MPI_Barrier(alltoall_comm_inter);
}
void glexcoll_Leader_baseda2a(int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	void *sendbuf = Get_sendbuf(bufmh, intra_rank);
	void *recvbuf = Get_recvbuf(bufmh, intra_rank);

	// puts("check 710");
	MPI_Barrier(MPI_COMM_WORLD);
	extern int Intra_Gather_Scatter_ON;
	// if (global_rank == 0)
	// 	puts("129");
	extern int num_of_ongoing_msg;
	//
	int block_size = sendsize * ppn;
	int shift = ppn * sendsize * alltoall_inter_rank;
	if (ppn != 1)
	{
		TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(shift, sendsize, bufmh);
	}
	else
	{
		memcpy(recvbuf + alltoall_inter_rank * sendsize, sendbuf + alltoall_inter_rank * sendsize, sendsize);
	}

	int message_length = block_size * ppn;
	int lenm = (1 << 27);
	int slice_count = message_length / lenm;
	if (message_length % lenm != 0)
		slice_count += 1;

	int total_length = block_size * ppn;
	int max_length = (1 << 27);
	int round_c = total_length / max_length;
	if (total_length % max_length != 0)
		round_c += 1;

	int inter_node_size = block_size * ppn;
	int steps_per_rounds = (BUFS_block_size) / inter_node_size;
	for (int h = 1; h < alltoall_inter_procn; h += steps_per_rounds)
	{

		int BufScount = 0;
		int start = h;
		int end = mmin(alltoall_inter_procn, start + steps_per_rounds);
		// if(alltoall_global_rank == 0) printf("start = %d end = %d round_c = %d\n",start,end,round_c);
		{
			char *p = bufS;
			for (int shift = start; shift < end; ++shift)
			{
				int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
				int leader = 0; // get_ith_leader(shift % leaderN);
				if (alltoall_intra_rank == leader)
				{
					glexcoll_intra_rma_gather_NUMA_Aware_with_mh(p, block_id * block_size, block_size, bufmh);
					p += block_size * ppn;
				}
			}
		}

		char *p = bufS;

		for (int shift = start; shift < end; ++shift)
		{
			int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
			int leader = 0; /// get_ith_leader(shift % leaderN);
			if (alltoall_intra_rank == leader)
			{
				int target = (alltoall_inter_rank + shift) % alltoall_inter_procn;
				int send_start = 0;
				for (int j = 0; j < round_c; j++)
				{
					int slice_length = mmin(max_length, total_length - send_start);
					// printf("rank = %d check 758\n",alltoall_global_rank);
					rdma_req.rmt_ep_addr.v = alltoall_ep_addrs[target].v; //_GLEXCOLL.ep_addrs[target].v;
					// printf("rank = %d check 760\n", alltoall_global_rank);
					rdma_req.local_mh.v = send_mh.v;
					rdma_req.local_offset = send_start + (p - bufS);
					rdma_req.len = slice_length;
					rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
					rdma_req.rmt_offset = send_start + BufScount * ppn * block_size;
					rdma_req.type = GLEX_RDMA_TYPE_PUT;
					rdma_req.rmt_evt.cookie[0] = 1 + alltoall_inter_rank;
					rdma_req.rmt_evt.cookie[1] = 1 + BufScount;
					rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
					if (j == round_c - 1)
					{
						if (BufScount % num_of_ongoing_msg == 1)
							rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
						else
							rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
					}
					else
						rdma_req.flag = NULL;
					rdma_req.next = NULL;
					int ret;
					while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
					{
					}
					if (ret != GLEX_SUCCESS)
					{
						if (ret == GLEX_INVALID_PARAM)
							printf("%d, _rdma() 非法参数", global_rank);
						printf("_rdma(), return: %d\n", ret);
						exit(1);
					}
					send_start += slice_length;
				}
				BufScount++;
				p += block_size * ppn;
			}
		}
		{
			//leaderèŠ‚ç‚¹å¿…é¡»ç­‰å¾…æ‰€æœ‰æ¶ˆæ¯åˆ°è¾
			for (int c = 0; c < BufScount; ++c)
			{
				while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
					;

				int source_inter_rank = event->cookie[0] - 1;
				int buf_recv_shift = event->cookie[1] - 1;
				//leaderèŠ‚ç‚¹è´Ÿè´£å°†æ¶ˆæ¯è½¬ç½
				int *source = bufR + buf_recv_shift * block_size * ppn;
				extern int MATRIX_Tanfform_ON;
				if (MATRIX_Tanfform_ON)
				{
					MATRIX_transform_local(source, ppn, sendsize);
				}
				//è½¬ç½®ä¹‹åŽï¼Œå°†æ¶ˆæ¯Scatteråˆ°å¯¹åº”å¤„
				// printf("%d recv from %d buf_recv_shift=%d\n",alltoall_inter_rank,source_inter_rank,buf_recv_shift);
				{
					glexcoll_intra_rma_scatter_NUMA_Aware_with_mh(source, source_inter_rank * block_size, block_size, bufmh);
				}
			}
		}
		if (BufScount > 0)
			glex_discard_probed_event(_GLEXCOLL.ep);
		if (alltoall_intra_rank == 0)
			MPI_Barrier(alltoall_comm_inter);
	}

	intra_memory_barrier_alltoall();
}

void rdma_send_recv_BUF(int start, int end, int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	int block_size = sendsize * ppn;
	int inter_node_sendsize = block_size * ppn;
	int num_of_ongoing_msg = 2;
	int max_length = (1 << 27);
	int round_c = inter_node_sendsize / max_length;
	if (inter_node_sendsize % max_length != 0)
		round_c += 1;
	int BufScount = 0;
	char *p = bufS;
	for (int shift = start; shift < end; shift++)
	{
		int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
		int leader = get_ith_leader(shift % leaderN);
		if (alltoall_intra_rank == leader)
		{
			int target = (alltoall_inter_rank + shift) % alltoall_inter_procn;
			int send_start = 0;
			for (int j = 0; j < round_c; j++)
			{
				int slice_length = mmin(max_length, inter_node_sendsize - send_start);
				// printf("rank = %d check 758\n",alltoall_global_rank);
				rdma_req.rmt_ep_addr.v = alltoall_ep_addrs[target].v; //_GLEXCOLL.ep_addrs[target].v;
				// printf("rank = %d check 760\n", alltoall_global_rank);
				rdma_req.local_mh.v = send_mh.v;
				rdma_req.local_offset = send_start + (p - bufS);
				rdma_req.len = slice_length;
				rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
				rdma_req.rmt_offset = send_start + BufScount * inter_node_sendsize;
				rdma_req.type = GLEX_RDMA_TYPE_PUT;
				rdma_req.rmt_evt.cookie[0] = 1 + alltoall_inter_rank;
				rdma_req.rmt_evt.cookie[1] = 1 + BufScount;
				rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
				if (j == round_c - 1)
				{
					if (BufScount % num_of_ongoing_msg == 1)
						rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
					else
						rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
				}
				else
					rdma_req.flag = NULL;
				rdma_req.next = NULL;
				int ret;
				while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
				{
				}
				if (ret != GLEX_SUCCESS)
				{
					if (ret == GLEX_INVALID_PARAM)
						printf("%d, _rdma() 非法参数", global_rank);
					printf("_rdma(), return: %d\n", ret);
					exit(1);
				}
				send_start += slice_length;
			}
			BufScount++;
			p += inter_node_sendsize;
		}
	}

	{
		for (int c = 0; c < BufScount; ++c)
		{
			while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
				;

			int source_inter_rank = event->cookie[0] - 1;
			int buf_recv_shift = event->cookie[1] - 1;
			//leaderèŠ‚ç‚¹è´Ÿè´£å°†æ¶ˆæ¯è½¬ç½
			int *source = bufR + buf_recv_shift * inter_node_sendsize;
			extern int MATRIX_Tanfform_ON;
			if (MATRIX_Tanfform_ON)
			{
				MATRIX_transform_local(source, ppn, sendsize);
			}
			// if(inter_rank == 0)
			//  	printf("%d recv from %d buf_recv_shift=%d\n",alltoall_inter_rank,source_inter_rank,buf_recv_shift);
			{
				glexcoll_intra_rma_scatter_NUMA_Aware_with_mh(source, source_inter_rank * block_size, block_size, bufmh);
			}
		}
	}
	if (am_i_leader())
	{
		glex_discard_probed_event(_GLEXCOLL.ep);
		// if(inter_rank == 0)
		// 	puts("1368");
		MPI_Barrier(alltoall_comm_inter);

		// if(inter_rank == 0)
		// 	puts("1370");
	}
}

void glexcoll_NMPML_a2a(int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)

{
	void *sendbuf = Get_sendbuf(bufmh, intra_rank);
	void *recvbuf = Get_recvbuf(bufmh, intra_rank);

	// puts("check 710");
	// MPI_Barrier(MPI_COMM_WORLD);
	extern int Intra_Gather_Scatter_ON;
	// puts("129");
	extern int num_of_ongoing_msg;

	int block_size = sendsize * ppn;
	int shift = ppn * sendsize * alltoall_inter_rank;
	if (ppn != 1)
	{
		TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(shift, sendsize, bufmh);
	}
	else
	{
		memcpy(recvbuf + alltoall_inter_rank * sendsize, sendbuf + alltoall_inter_rank * sendsize, sendsize);
	}

	int inter_node_size = block_size * ppn;
	int steps_per_rounds = leaderN * ((BUFS_block_size) / inter_node_size);

	for (unsigned int h = 1; h < alltoall_inter_procn; h += steps_per_rounds)
	{

		int start = h;
		int end = mmin(alltoall_inter_procn, start + steps_per_rounds);
		// if(global_rank == 0)
		// 	printf("start= %d end = %d\n",start,end);

		{
			char *p = bufS;
			for (int shift = start; shift < end; ++shift)
			{
				int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
				int leader = get_ith_leader(shift % leaderN);
				if (alltoall_intra_rank == leader)
				{
					if (p - bufS >= BUFS_block_size)
					{
						printf("rdma发送缓冲区溢出,global rank = %d shift=%d leader =%d\n", global_rank, shift, leader);
						// exit(0);
					}
					glexcoll_intra_rma_gather_NUMA_Aware_with_mh(p, block_id * block_size, block_size, bufmh);
					p += inter_node_size;
					// if(global_rank == 8){
					// 	printf("shift = %d,sz=%d\n",shift,inter_node_size);
					// }
				}
			}
		}

		rdma_send_recv_BUF(start, end, sendsize, bufmh);
	}

	intra_memory_barrier_alltoall();
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(global_rank == 0)
	// 	puts("NMPML finish");
}

void glexcoll_multileader_pipeline_rdma_Topology_aware_alltoall_with_mh(int sendsize,
																		struct GLEXCOLL_a2a_bufmh *bufmh)
{
	void *sendbuf = Get_sendbuf(bufmh, intra_rank);
	void *recvbuf = Get_recvbuf(bufmh, intra_rank);

	// MPI_Barrier(MPI_COMM_WORLD);
	// if(alltoall_global_rank == 0) puts("check 849");
	extern int Intra_Gather_Scatter_ON;
	// puts("129");
	extern int num_of_ongoing_msg;
	//ç¬¬ä¸€æ­¥å…ˆè¿›è¡ŒèŠ‚ç‚¹å†…alltoallè½¬ç½®
	int block_size = sendsize * ppn;
	int shift = ppn * sendsize * alltoall_inter_rank;
	if (ppn != 1)
	{
		TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(shift, sendsize, bufmh);
	}
	else
	{
		memcpy(recvbuf + alltoall_inter_rank * sendsize, sendbuf + alltoall_inter_rank * sendsize, sendsize);
	}
	// MPI_Barrier(MPI_COMM_WORLD);    if(alltoall_global_rank == 0) puts("check 721");

	// printf("sendbuf[myself]=%d, recvbuf[myself=%d]=%d\n", ((int *)sendbuf)[alltoall_inter_rank], alltoall_inter_rank, ((int *)recvbuf)[alltoall_inter_rank]);
	// TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(shift, sendsize);
	//ç¬¬äºŒæ­¥è¿›è¡ŒèŠ‚ç‚¹é—´alltoallè½¬ç½®
	//é¦–å…ˆè¿›è¡ŒGatheræ“ä½œï¼Œå°†æ•°æ®ç”±å¤šä¸ªLeaderå„è‡ª Gatheråˆ°è‡ªå·±çš„BufS

	int BufScount = 0;
	{
		char *p = bufS;
		for (int index = 1; index < alltoall_inter_procn; ++index)
		{
			// MPI_Barrier(MPI_COMM_WORLD);    if(alltoall_global_rank == 0) puts("check 869");
			shift = (alltoall_inter_procn + alltoall_send_order[index] - alltoall_inter_rank) % alltoall_inter_procn;
			// MPI_Barrier(MPI_COMM_WORLD);    if(alltoall_global_rank == 0) puts("check 870");

			int block_id = (alltoall_inter_rank + shift) % alltoall_inter_procn;
			//leaderNæ˜¯ç†è®ºå€
			int ActualLeaderN = mmax(mmin(ppn / 4, 4), 1);
			int leader = (shift % ActualLeaderN) << 2;
			// printf("leader = %d\n", leader);
			// if(intra_rank == 0 && alltoall_inter_rank == 7)
			//     printf("alltoall_inter_rank=%d's shift = %d alltoall_send_order[index=%d]=%d\n", alltoall_inter_rank, shift,index,alltoall_send_order[index]);
			if (alltoall_intra_rank == leader)
			{
				// printf("leader = %d\n", leader);
				glexcoll_intra_rma_gather_with_mh(p, block_id * block_size, block_size, bufmh);
				{
					// æ‰“å°gatherè¾“å‡ºæ•°æ®
					// if(intra_rank == 0 && alltoall_inter_rank == 7)
					// {
					//     // puts("\n-------------------------------");
					//     int i = 0;
					//     for(int a = 0;a<ppn;a++)
					//         for(int s=0;s<ppn;s++)
					//     {
					//         int source = alltoall_global_rank - (alltoall_global_rank % ppn) + a;
					//         int target = ((alltoall_inter_rank + shift) % alltoall_inter_procn)*ppn + s;
					//         int val = source*100000 + target;
					//         // if(((((int *)p)[i]) % 10000000) != val)
					//         {
					//             printf("ERROR : source=%d target = %d gather result: global_rank=%d %d\n",source,target, global_rank, ((int *)p)[i]);
					//         }
					//         i++;
					//     }
					// }
				}
				// printf("%d finish glexcoll_intra_rma_gather_with_mh\n", alltoall_global_rank, shift);
				// MPI_Barrier(alltoall_comm_inter);
				//æ•°æ®gatherå®ŒæˆåŽç«‹é©¬å‘é€å‡ºåŽ
				{
					int target = alltoall_send_order[index]; //(alltoall_inter_rank + shift) % alltoall_inter_procn;
					// printf("rank = %d target inter_rank = %d\n",alltoall_global_rank,target);
					rdma_req.rmt_ep_addr.v = alltoall_ep_addrs[target].v; //_GLEXCOLL.ep_addrs[target].v;
					// printf("rank = %d check 760\n", alltoall_global_rank);
					rdma_req.local_mh.v = send_mh.v;
					rdma_req.local_offset = (p - bufS);
					rdma_req.len = block_size * ppn;
					rdma_req.rmt_mh.v = _GLEXCOLL.rmt_mhs[target].v;
					rdma_req.rmt_offset = (alltoall_inter_rank / ActualLeaderN) * ppn * block_size;
					rdma_req.type = GLEX_RDMA_TYPE_PUT;
					rdma_req.rmt_evt.cookie[0] = 1 + alltoall_inter_rank;
					rdma_req.rmt_evt.cookie[1] = 1 + alltoall_inter_rank / ActualLeaderN;
					rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
					if (BufScount % 2 == 1)
						rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
					else
						rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
					// rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;

					rdma_req.next = NULL;
					int ret;
					while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
					{
					}
					if (ret != GLEX_SUCCESS)
					{
						if (ret == GLEX_INVALID_PARAM)
							printf("%d, _rdma() 非法参数", global_rank);
						printf("_rdma(), return: %d\n", ret);
						exit(1);
					}
					BufScount++;
				}
				//å‘é€æŒ‡é’ˆæ›´æ–
				p += block_size * ppn;
			}
		}
	}
	//  MPI_Barrier(MPI_COMM_WORLD);    if(alltoall_global_rank == 0) puts("check 788");

	{
		//leaderèŠ‚ç‚¹å¿…é¡»ç­‰å¾…æ‰€æœ‰æ¶ˆæ¯åˆ°è¾
		for (int c = 0; c < BufScount; ++c)
		{
			// printf("%d's BufScount = %d\n",alltoall_global_rank,BufScount);
			while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
				;

			int source_inter_rank = event->cookie[0] - 1;
			int buf_recv_shift = event->cookie[1] - 1;
			//leaderèŠ‚ç‚¹è´Ÿè´£å°†æ¶ˆæ¯è½¬ç½
			int *source = bufR + buf_recv_shift * block_size * ppn;
			extern int MATRIX_Tanfform_ON;
			if (MATRIX_Tanfform_ON)
			{
				MATRIX_transform_local(source, ppn, sendsize);
			}
			//è½¬ç½®ä¹‹åŽï¼Œå°†æ¶ˆæ¯Scatteråˆ°å¯¹åº”å¤„
			// printf("%d recv from %d buf_recv_shift=%d data=%d\n",alltoall_global_rank,source_inter_rank,buf_recv_shift,*source);
			{
				glexcoll_intra_rma_scatter_with_mh(source, source_inter_rank * block_size, block_size, bufmh);
			}
			//  puts("finish glexcoll_intra_rma_scatter_with_mh");
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// puts("check 942");
	if (BufScount > 0)
		glex_discard_probed_event(_GLEXCOLL.ep);
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(alltoall_global_rank == 0)
	//     puts("check 946");

	intra_memory_barrier_alltoall();
	if (alltoall_intra_rank == 0)
		MPI_Barrier(alltoall_comm_inter);
	// printf("sendbuf[myself]=%d, recvbuf[myself=%d]=%d\n", ((int *)sendbuf)[alltoall_inter_rank], alltoall_inter_rank, ((int *)recvbuf)[alltoall_inter_rank]);
}

//将alltoall切分为多次传输

#pragma optimize("", off)
void load_to_cache(char *addr, int sz)
{
	static char v;
	for (int i = 0; i < sz; i++)
	{
		v = addr[i];
		addr[i] = v;
	}
}
#pragma optimize("", on)

void glexcoll_Split_based_alltoall_to_multiple_with_mh_ppn1(int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	int tmp = sendsize / 16;
	int slice_sz = mmax(8192, tmp);
	slice_sz = 8192;
	int slice_count = sendsize / slice_sz;
	if (slice_count < 1)
		slice_count = 1;
	int slice_shift = 0;

	int remain_sz = sendsize;
	for (int h = 0; h < slice_count; h++, slice_shift += slice_sz)
	{
		int realsize = mmin(remain_sz, slice_sz);

		int shift = alltoall_inter_rank * sendsize;
		int BufScount = 0;
		for (int i = 0; i < alltoall_inter_procn; i++)
		{
			int target;
			target = (alltoall_inter_rank + i) % alltoall_inter_procn;
			if (target == alltoall_inter_rank)
				memcpy(bufmh->recvbuf + shift, bufmh->sendbuf + shift, sendsize);
			else
			{
				int local_offset = target * sendsize + slice_shift;
				// load_to_cache(bufmh->sendbuf + local_offset, realsize);
				int rmt_offset = alltoall_inter_rank * sendsize + slice_shift;
				// load_to_cache(bufmh->recvbuf + rmt_offset, realsize);

				rdma_req.rmt_ep_addr.v = (bufmh->alltoall_ep_addrs[target]).v;
				rdma_req.local_mh.v = (bufmh->sendmh).v;
				rdma_req.local_offset = local_offset;
				rdma_req.len = realsize;
				rdma_req.rmt_mh.v = (bufmh->alltoall_mem_handle_vec[target]).v;
				rdma_req.rmt_offset = rmt_offset;
				rdma_req.type = GLEX_RDMA_TYPE_PUT;
				rdma_req.rmt_evt.cookie[0] = 1 + 1;
				rdma_req.rmt_evt.cookie[1] = 1 + 1;
				rdma_req.local_evt.cookie[0] = 1 + 9;
				rdma_req.local_evt.cookie[1] = 1 + 9;
				rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
				rdma_req.flag = NULL;
				// if (BufScount % 2 == 1)
				//     rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
				// else
				//     rdma_req.flag =  GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
				rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
				rdma_req.next = NULL;
				int ret;
				while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
				{
				}
				if (ret != GLEX_SUCCESS)
				{
					if (ret == GLEX_INVALID_PARAM)
						printf("%d, _rdma() 非法参数", global_rank);
					printf("_rdma(), return: %d\n", ret);
					exit(1);
				}
				BufScount++;
			}
		}
		for (int i = 1; i < alltoall_inter_procn; i++)
		{
			{
				while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
					;
				// printf("get a event id = %d my-rank =%d \n",event->cookie[1], inter_rank);

				extern void PJT_discard_event(glex_ep_handle_t ep);
				PJT_discard_event(_GLEXCOLL.ep);
			}
		}
	}
	// for (int i = 1; i < alltoall_inter_procn; i++)
	//     for (int h = 0; h < slice_count; h++)
	//  MPI_Barrier(MPI_COMM_WORLD);
}
//将一个大消息a2a切分为多个。然后降低event的使用，增加并发度。
void glexcoll_Split_based_alltoall_to_multiple_with_mh_ppn1_event_reduce(int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	int tmp = sendsize / 16;
	int slice_sz = mmax(8192, tmp);
	slice_sz = 8192;
	int slice_count = sendsize / slice_sz;
	if (slice_count < 1)
		slice_count = 1;
	int slice_shift = 0;

	int remain_sz = sendsize;
	for (int h = 0; h < slice_count; h++, slice_shift += slice_sz)
	{
		int realsize = mmin(remain_sz, slice_sz);

		int shift = alltoall_inter_rank * sendsize;
		int BufScount = 0;
		for (int i = 0; i < alltoall_inter_procn; i++)
		{
			int target;
			target = (alltoall_inter_rank + i) % alltoall_inter_procn;
			if (target == alltoall_inter_rank)
				memcpy(bufmh->recvbuf + shift, bufmh->sendbuf + shift, sendsize);
			else
			{
				int local_offset = target * sendsize + slice_shift;
				int rmt_offset = alltoall_inter_rank * sendsize + slice_shift;
				rdma_req.rmt_ep_addr.v = (bufmh->alltoall_ep_addrs[target]).v;
				rdma_req.local_mh.v = (bufmh->sendmh).v;
				rdma_req.local_offset = local_offset;
				rdma_req.len = realsize;
				rdma_req.rmt_mh.v = (bufmh->alltoall_mem_handle_vec[target]).v;
				rdma_req.rmt_offset = rmt_offset;
				rdma_req.type = GLEX_RDMA_TYPE_PUT;
				rdma_req.rmt_evt.cookie[0] = 1 + 1;
				rdma_req.rmt_evt.cookie[1] = 1 + 1;
				rdma_req.local_evt.cookie[0] = 1 + 9;
				rdma_req.local_evt.cookie[1] = 1 + 9;
				rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
				rdma_req.flag = NULL;
				// if (BufScount % 2 == 1)
				//     rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
				// else
				//     rdma_req.flag =  GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
				rdma_req.flag = NULL;
				if (i == alltoall_inter_procn - 1)
					rdma_req.flag = GLEX_FLAG_LOCAL_EVT | GLEX_FLAG_FENCE;
				rdma_req.next = NULL;
				int ret;
				while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
				{
				}
				if (ret != GLEX_SUCCESS)
				{
					if (ret == GLEX_INVALID_PARAM)
						printf("%d, _rdma() 非法参数", global_rank);
					printf("_rdma(), return: %d\n", ret);
					exit(1);
				}
				BufScount++;
			}
		}
		// for (int i = 1; i < alltoall_inter_procn; i++)
		for (int i = 0; i < 1; i++)
		{
			{
				while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
					;
				extern void PJT_discard_event(glex_ep_handle_t ep);
				PJT_discard_event(_GLEXCOLL.ep);
			}
		}
	}
}

//将一个大消息a2a切分为多个。然后降低event的使用，增加并发度。
int get_my_target_at_round(int rank, int procn, int round, int i)
{
	int get_position(int procn, int round, int i)
	{
		if (i < procn - round)
			return round + i;
		else
		{
			return procn - i - 1;
		}
	}
	return get_position(procn, round, (rank + i) % procn);
}
int alltoall_reorder_round = 0;
void glexcoll_Split_based_alltoall_to_multiple_with_mh_ppn1_event_reduce_reorder(int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	int tmp = sendsize / 16;
	int slice_sz = mmax(8192, tmp);
	slice_sz = 8192;
	int slice_count = sendsize / slice_sz;
	if (slice_count < 1)
		slice_count = 1;
	int slice_shift = 0;

	int remain_sz = sendsize;
	for (int h = 0; h < slice_count; h++, slice_shift += slice_sz)
	{
		int realsize = mmin(remain_sz, slice_sz);

		int shift = alltoall_inter_rank * sendsize;
		int BufScount = 0;
		for (int i = 0; i < alltoall_inter_procn; i++)
		{
			int a = i / 2;
			int b = i & 0x1;

			int target = get_my_target_at_round(alltoall_inter_rank, alltoall_inter_procn, alltoall_reorder_round, i);
			// if(b == 0)
			//     target = (alltoall_inter_rank + a) % alltoall_inter_procn;
			// else
			//     target = (alltoall_inter_rank + alltoall_inter_procn - 1 -a) % alltoall_inter_procn;
			if (target == alltoall_inter_rank)
				memcpy(bufmh->recvbuf + shift, bufmh->sendbuf + shift, sendsize);
			else
			{
				int local_offset = target * sendsize + slice_shift;
				int rmt_offset = alltoall_inter_rank * sendsize + slice_shift;
				rdma_req.rmt_ep_addr.v = (bufmh->alltoall_ep_addrs[target]).v;
				rdma_req.local_mh.v = (bufmh->sendmh).v;
				rdma_req.local_offset = local_offset;
				rdma_req.len = realsize;
				rdma_req.rmt_mh.v = (bufmh->alltoall_mem_handle_vec[target]).v;
				rdma_req.rmt_offset = rmt_offset;
				rdma_req.type = GLEX_RDMA_TYPE_PUT;
				rdma_req.rmt_evt.cookie[0] = 1 + 1;
				rdma_req.rmt_evt.cookie[1] = 1 + 1;
				rdma_req.local_evt.cookie[0] = 1 + 9;
				rdma_req.local_evt.cookie[1] = 1 + 9;
				rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
				rdma_req.flag = NULL;
				// if (BufScount % 2 == 1)
				//     rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
				// else
				//     rdma_req.flag =  GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
				rdma_req.flag = NULL;
				if (i == alltoall_inter_procn - 1)
					rdma_req.flag = GLEX_FLAG_LOCAL_EVT | GLEX_FLAG_FENCE;
				rdma_req.next = NULL;
				int ret;
				while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
				{
				}
				if (ret != GLEX_SUCCESS)
				{
					if (ret == GLEX_INVALID_PARAM)
						printf("%d, _rdma() 非法参数", global_rank);
					printf("_rdma(), return: %d\n", ret);
					exit(1);
				}
				BufScount++;
			}
		}
		// for (int i = 1; i < alltoall_inter_procn; i++)
		for (int i = 0; i < 1; i++)
		{
			{
				while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
					;
				extern void PJT_discard_event(glex_ep_handle_t ep);
				PJT_discard_event(_GLEXCOLL.ep);
			}
		}
	}
}

void glexcoll_event_avoid_alltoall_with_mh_ppn1(int sendsize,
												struct GLEXCOLL_a2a_bufmh *bufmh)
{
	int shift = alltoall_inter_rank * sendsize;
	int BufScount = 0;
	int Bias = 0;
	if (alltoall_inter_rank % 2 == 0)
		Bias = 0;
	else
		Bias = 0; //alltoall_inter_procn / 2;
	for (int i = 0; i < alltoall_inter_procn; i++)
	{
		int target;
		target = (alltoall_inter_rank + Bias + i) % alltoall_inter_procn; //alltoall_send_order[i];
		//target = alltoall_send_order[i];
		// if(inter_procn <=64)
		//     target = (alltoall_inter_rank + Bias + i) % alltoall_inter_procn;//alltoall_send_order[i];
		// else
		//     target = alltoall_send_order[i];
		if (target == alltoall_inter_rank)
			memcpy(bufmh->recvbuf + shift, bufmh->sendbuf + shift, sendsize);
		else
		{
			rdma_req.rmt_ep_addr.v = (bufmh->alltoall_ep_addrs[target]).v; //_GLEXCOLL.ep_addrs[target].v;
			// printf("rank = %d target = %d check 760\n", alltoall_global_rank,target);
			rdma_req.local_mh.v = (bufmh->sendmh).v;
			rdma_req.local_offset = target * sendsize;
			rdma_req.len = sendsize;
			rdma_req.rmt_mh.v = (bufmh->alltoall_mem_handle_vec[target]).v;
			rdma_req.rmt_offset = alltoall_inter_rank * sendsize;
			rdma_req.type = GLEX_RDMA_TYPE_PUT;
			rdma_req.rmt_evt.cookie[0] = 1 + 1;
			rdma_req.rmt_evt.cookie[1] = 1 + 1;
			rdma_req.local_evt.cookie[0] = 1 + 6;
			rdma_req.local_evt.cookie[1] = 1 + 6;
			rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;

			rdma_req.flag = NULL;
			if (BufScount % 2 == 1)
				rdma_req.flag = GLEX_FLAG_FENCE;

			if (i == alltoall_inter_procn - 1)
				rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_LOCAL_EVT | GLEX_FLAG_FENCE;
			else
				// if (i == alltoall_inter_procn - 1)
				//     rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_LOCAL_EVT | GLEX_FLAG_FENCE;
				// else{
				//     rdma_req.flag = NULL;
				// }
				rdma_req.next = NULL;
			int ret;
			while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
			{
			}
			if (ret != GLEX_SUCCESS)
			{
				if (ret == GLEX_INVALID_PARAM)
					printf("%d, _rdma() 非法参数", global_rank);
				printf("_rdma(), return: %d\n", ret);
				exit(1);
			}
			BufScount++;
		}
	}

	// for(int i = 1;i< alltoall_inter_procn;i++)
	for (int j = 0; j < 2; j++)
	{
		{
			// printf("%d's BufScount = %d\n",alltoall_global_rank,BufScount);
			while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
				;
			// printf("get a event id = %d my-rank =%d \n",event->cookie[1], inter_rank);
			extern void PJT_discard_event(glex_ep_handle_t ep);
			PJT_discard_event(_GLEXCOLL.ep);
			//  puts("finish glexcoll_intra_rma_scatter_with_mh");
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
void glexcoll_RDMA_alltoall_with_mh_ppn1(int sendsize, struct GLEXCOLL_a2a_bufmh *bufmh)
{
	extern int num_of_ongoing_msg;
	int shift = alltoall_inter_rank * sendsize;
	int BufScount = 0;
	for (int i = 0; i < alltoall_inter_procn; i++)
	{
		if (i == 0)
			memcpy(bufmh->recvbuf + shift, bufmh->sendbuf + shift, sendsize);
		else
		{
			int target;
			extern int Topology_aware_up;
			target = (alltoall_inter_rank + i) % alltoall_inter_procn;
			rdma_req.rmt_ep_addr.v = (bufmh->alltoall_ep_addrs[target]).v; //_GLEXCOLL.ep_addrs[target].v;
			// printf("rank = %d target = %d check 760\n", alltoall_global_rank,target);
			rdma_req.local_mh.v = (bufmh->sendmh).v;
			rdma_req.local_offset = target * sendsize;
			rdma_req.len = sendsize;
			rdma_req.rmt_mh.v = (bufmh->alltoall_mem_handle_vec[target]).v;
			rdma_req.rmt_offset = alltoall_inter_rank * sendsize;
			rdma_req.type = GLEX_RDMA_TYPE_PUT;
			rdma_req.rmt_evt.cookie[0] = 1 + 2;
			rdma_req.rmt_evt.cookie[1] = 1 + 2;
			rdma_req.rmt_key = _GLEXCOLL.ep_attr.key;
			if (BufScount % num_of_ongoing_msg == 0)
				rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
			else
				rdma_req.flag = GLEX_FLAG_REMOTE_EVT;
			// rdma_req.flag = GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
			rdma_req.next = NULL;
			int ret;
			while ((ret = glex_rdma(_GLEXCOLL.ep, &rdma_req, &bad_rdma_req)) == GLEX_BUSY)
			{
			}
			if (ret != GLEX_SUCCESS)
			{
				if (ret == GLEX_INVALID_PARAM)
					printf("%d, _rdma() 非法参数", global_rank);
				printf("_rdma(), return: %d\n", ret);
				exit(1);
			}
			BufScount++;
		}
	}
	for (int i = 1; i < alltoall_inter_procn; i++)
	{
		{
			// printf("%d's BufScount = %d\n",alltoall_global_rank,BufScount);
			while (glex_probe_next_event(_GLEXCOLL.ep, &event) == GLEX_NO_EVENT)
				;
			extern void PJT_discard_event(glex_ep_handle_t ep);
			PJT_discard_event(_GLEXCOLL.ep);
			//  puts("finish glexcoll_intra_rma_scatter_with_mh");
		}
	}
}
extern void DIRECT_alltoall(void *sendbuf,
							int sendsize,
							void *recvbuf,
							int recvsize,
							MPI_Comm comm,
							MPI_Datatype type);

int Topology_aware_up = 0;
void DIRECT_Topology_aware_alltoall(void *sendbuf,
									int sendsize,
									void *recvbuf,
									int recvsize,
									MPI_Comm comm,
									MPI_Datatype type)
{
	int alltoall_rank;
	int alltoall_procn;
	MPI_Comm_rank(comm, &alltoall_rank);
	MPI_Comm_size(comm, &alltoall_procn);
	static MPI_Request *reqSV;	//[128];
	static MPI_Request *reqRV;	//[128];
	static MPI_Status *statusV; //[128];
	if (reqSV == 0)
		reqSV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoall_procn);
	if (reqRV == 0)
		reqRV = (MPI_Request *)malloc(sizeof(MPI_Request) * alltoall_procn);
	if (statusV == 0)
		statusV = (MPI_Status *)malloc(sizeof(MPI_Status) * alltoall_procn);
	//MPI_Alltoall(sendbuf,sendsize,MPI_DOUBLE,recvbuf,recvsize,MPI_DOUBLE,comm);
	// MPI_Request req_sends[global_procn];
	// MPI_Request req_recvs[global_procn];
	// MPI_Status status[global_procn];
	// if(global_rank ==1)
	//     printf("sendbuf[0] = %f,sendbuf[1] = %f\n",((double *)sendbuf)[0],((double *)sendbuf)[1]);
	int s = 0;
	int SMAX = alltoall_procn;
	for (int shiftb = 0; shiftb < alltoall_procn; shiftb += SMAX)
	{
		s = 0;
		for (s = 0; s < SMAX && shiftb + s < alltoall_procn; s++)
		// if(shift >=0 && shift < 5)
		{
			int shift = shiftb + s;
			// printf("shift = %d\n",shift);
			{
				int target, source;
				if (Topology_aware_up == 1)
				{
					target = alltoall_send_order[shift];
					source = alltoall_recv_order[shift];
					//  if(shift == 1)
					//  printf("%d's target = %d\n",alltoall_inter_rank,target);
				}
				else
				{
					target = (alltoall_rank + shift) % alltoall_procn;
					source = (alltoall_procn + alltoall_rank - shift) % alltoall_procn;
				}
				// if(alltoall_rank == 0 || alltoall_rank == 4)
				// printf("%d's target = %d\n",alltoall_inter_rank,target);
				MPI_Isend(sendbuf + sendsize * target, sendsize, MPI_CHAR,
						  target, alltoall_rank, comm, &(reqSV[s]));
				MPI_Irecv(recvbuf + recvsize * source, sendsize, MPI_CHAR,
						  source,
						  source, comm, &(reqRV[s]));
				// MPI_Wait(&(reqSV[shift]), statusV);
				// MPI_Wait(&(reqRV[shift]), statusV);
			}
		}
		MPI_Waitall(s, reqSV, statusV);
		MPI_Waitall(s, reqRV, statusV);
		// MPI_Barrier(comm);
		// if(alltoall_rank == 0) puts("---------------------------------------------------------------------------------------");
		// MPI_Barrier(comm);
		// else{
		//         int target,source;
		//         if(Topology_aware_up == 1){
		//              source = target = alltoall_send_order[shift];
		//             //  if(shift == 1)
		//             //  printf("%d's target = %d\n",alltoall_inter_rank,target);
		//         }
		//         else                       {
		//             target = (alltoall_rank + shift) % alltoall_procn;
		//             source = (alltoall_procn + alltoall_rank - shift) % alltoall_procn;
		//         }
		//         // if(alltoall_rank == 0 || alltoall_rank == 4)
		//         // printf("%d's target = %d\n",alltoall_inter_rank,target);
		//         MPI_Isend(sendbuf + sendsize * target, sendsize, MPI_CHAR,
		//                   target, alltoall_rank, comm, &(reqSV[s]));
		//         MPI_Irecv(recvbuf + recvsize * source, sendsize, MPI_CHAR,
		//                   source,
		//                   source, comm, &(reqRV[s]));
		//         s++;
		// }
	}
	// int shift = 5;
	// int source,target;
	//             if(Topology_aware_up == 1){
	//                 shift = 1;
	//                  source = target = alltoall_send_order[shift];
	//                 //  if(shift == 1)
	//                 //  printf("%d's target = %d\n",alltoall_inter_rank,target);
	//             }
	//             else                       {
	//                 target = (alltoall_rank + shift) % alltoall_procn;
	//                 source = (alltoall_procn + alltoall_rank - shift) % alltoall_procn;
	//             }
	//             if(alltoall_rank == 0 || alltoall_rank == 6)
	//                 printf("%d's target = %d\n",alltoall_inter_rank,target);
	//             MPI_Isend(sendbuf + sendsize * target, sendsize, MPI_CHAR,
	//                       target, alltoall_rank, comm, &(reqSV[shift]));
	//             MPI_Irecv(recvbuf + recvsize * source, sendsize, MPI_CHAR,
	//                       source,
	//                       source, comm, &(reqRV[shift]));
	//             MPI_Wait(&(reqSV[shift]), statusV);
	//             MPI_Wait(&(reqRV[shift]), statusV);
}
void GLEXCOLL_Alltoall_new(void *sendbuf,
						   int sendsize,
						   void *recvbuf,
						   int recvsize,
						   MPI_Comm comm,
						   MPI_Datatype type)
{
	//å¯¹äºŽä¸€ä¸ªèŠ‚ç‚¹çš„æƒ…å†µ

	DIRECT_Topology_aware_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm, type);
	// if (alltoall_inter_procn == 1)
	// {
	//     TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer(0, sendsize);
	//     return 0;
	// }
	// else if (alltoall_global_procn <= 4)
	// {
	//     DIRECT_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm, type);
	//     return 0;
	// }
	// else
	// {
	//     {
	//         glexcoll_multileader_pipeline_rdma_alltoall(sendbuf, sendsize, recvbuf, recvsize, comm);
	//     }
	// }
}

void GLEXCOLL_Alltoall_pjt(struct GLEXCOLL_a2a_bufmh *bufmh, int size)
{

	// Topology_aware_up = 0;
	// Alltoall_algorithm = ONMPML;
	// if(alltoall_global_rank == 0)
	// {
	//     printf("glexcoll: chech 2135 size = %d\n",size);
	// }
	if (alltoall_inter_procn == 1)
	{
		TPDS17_Cache_oblivious_intra_node_NUMA_alltoall_on_buffer_new(0, size, bufmh);
		return 0;
	}
	else
	{
		if (ppn == 1)
		{
			switch (Topology_aware_up)
			{
			case 0:
				glexcoll_RDMA_alltoall_with_mh_ppn1(size, bufmh);
				break;
			case 1:
				glexcoll_Split_based_alltoall_to_multiple_with_mh_ppn1(size, bufmh);
				break;
			case 2:
				glexcoll_Split_based_alltoall_to_multiple_with_mh_ppn1_event_reduce(size, bufmh);
				break;
			case 3:
				glexcoll_Split_based_alltoall_to_multiple_with_mh_ppn1_event_reduce_reorder(size, bufmh);
				break;
			}
		}
		else
		{
			if (Topology_aware_up == 1)
			{
				glexcoll_multileader_pipeline_rdma_Topology_aware_alltoall_with_mh(size, bufmh);
			}
			else
				switch (Alltoall_algorithm)
				{
				case ONMPML:
					glexcoll_ONMPML_a2a(size, bufmh);
					break;
				case L_a2a:
					glexcoll_Leader_baseda2a(size, bufmh);
					break;
				case NMPML:
					glexcoll_NMPML_a2a(size, bufmh);
					break;
				case SONMPML:
					glexcoll_SONMPML_a2a_with_mh(size, bufmh);
				}
		}
	}
}
void GLEXCOLL_Alltoall_Finalize_new()
{
	extern void glexcoll_destroy_alltoall_shared_memory_buffer();
	glexcoll_destroy_alltoall_shared_memory_buffer();

#ifdef PJT_NEW_VERSION
	if (alltoall_inter_procn > 1 && am_i_leader())
#else
	if (alltoall_inter_procn > 1 && alltoall_intra_rank % 4 == 0 && (alltoall_intra_rank >> 2) < leaderN)
#endif
	{

		int ret = glex_deregister_mem(_GLEXCOLL.ep, _GLEXCOLL.local_mh);
		if (ret != GLEX_SUCCESS)
		{
			printf("_deregister_mem(), return: %d\n", ret);
			exit(1);
		}

		ret = glex_deregister_mem(_GLEXCOLL.ep, send_mh);
		if (ret != GLEX_SUCCESS)
		{
			printf("_deregister_mem(), return: %d\n", ret);
			exit(1);
		}
	}

	free(bufS);
	free(bufR);
}

/*
   int main()
   {
   int procn = 10;
   for(int round = 0;round < 10;round++)
   {

   for(int i = 0;i<10;i++)
   {
   printf("%d ",get_my_target_at_round(0,procn,round,i));
   }
   puts("");
   }
   return 0;
   }
   */
