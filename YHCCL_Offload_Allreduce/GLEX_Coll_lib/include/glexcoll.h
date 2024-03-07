#ifndef GLEXCOLL_H
#define GLEXCOLL_H
#include <mpi.h>
#include <pthread.h>
#include <semaphore.h>
#include "glex.h"

#define PERFORMANCE_DEBUG

#ifdef __cplusplus
#define BEGIN_C_DECLS \
	extern "C"        \
	{
#define END_C_DECLS }
#else /* !__cplusplus */
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif /* __cplusplus */

BEGIN_C_DECLS
struct Payload
{
	char lock;
	char flag;
	int size;
	union
	{
		double data_double[7];
		int data_int[7];
		float data_float[7];
	} T;
};

typedef union
{
	char V[64];
	struct Payload payload;
} SHM_FLIT;

#define EP_CREDIT_MAX 128
#define Event_CREDIT_MAX 128
struct _GLEXCOLL_Data
{
	uint32_t num_of_devices;
	struct glex_ep_attr ep_attr, ep_attr1;
	glex_device_handle_t dev;
	struct glex_device_attr dev_attr;
	glex_ep_handle_t ep;
	int ep_credit;
	int event_credit;
	glex_ep_handle_t ep1;
	glex_ep_addr_t *ep_addrs;
	glex_ep_addr_t *ep_addrs1;
	glex_mem_handle_t local_mh;
	glex_mem_handle_t *rmt_mhs;

	MPI_Comm Comm_Global;
	int global_procn;
	int global_rank;
	MPI_Comm Comm_intra;
	int intra_procn;
	int intra_rank;
	MPI_Comm Comm_inter;
	int inter_procn;
	int inter_rank;
	MPI_Comm Comm_numa_intra;
	int numa_intra_procn;
	int numa_intra_rank;
	MPI_Comm Comm_numa_inter;
	int numa_inter_procn;
	int numa_inter_rank;
	char host_name[MPI_MAX_PROCESSOR_NAME];
	int *ReduceId_to_InterRank_Map;
	struct GLEXCOLL_ALLREDUCE_PARA
	{
		void *sendbuf;
		void *recvbuf;
		MPI_Datatype datatype;
		enum glex_coll_reduce_data_bits bit_size;
		int count;
		MPI_Op op;
		int root; //allreduce_root
	} _allreduce;
	//非阻塞同步所使用的互斥锁和环境变量
	pthread_cond_t cond;
	pthread_mutex_t mtx;
	pthread_t comm_thread;
	sem_t PostREQ, FinishREQ;
	enum threadFLAG
	{
		RUNNING,
		EXIT,
		ALLREDUCE,
		REDUCE,
		BCAST
	} th_flag;
	//共享内存所需的变量
	char name[50];
	volatile char *SHM_buf;
	volatile char *SHM_0;
	volatile void *SHM_bufs[64];
};
enum node_type
{
	LEAF,
	MID,
	ROOT
};
struct TreeLevelInfo
{
	int parentID;
	glex_ep_addr_t parentAddr;		 //有效数据
	int brothersN;					 //有效数据
	int UpReduceRank;				 //仅在LEAF和MID点中有效，上方规约参与点中的rank，
	int DownReduceRank;				 //仅在ROOT和MID点中有效，					  //有效数据
	int childsN;					 //有效数据
	enum node_type type;			 //有效数据
	glex_ep_addr_t childAddrs[256];	 //有效数据
	glex_ep_addr_t childAddrs1[256]; //有效数据
	int childIds[256];
	enum node_type childTypes[256]; //最新有效数据，
};

enum INTER_COMM_MODE
{
	OFFLOAD,
	NON_OFFLOAD
};
extern int ppn;
extern MPI_Comm Comm_intra;
extern int intra_procn;
extern int intra_rank;
extern MPI_Comm Comm_inter;
extern int inter_procn;
extern int inter_rank;
extern int global_rank;
extern int global_procn;
extern struct _GLEXCOLL_Data _GLEXCOLL;
extern union glex_coll_data_loc data_loc[GLEX_COLL_REDUCE_MAX_DATA_UNITS];
extern struct glex_imm_mp_req GLEX_Coll_req, GLEX_Coll_req_backup;
extern char GLEXCOLL_databuf[1024];
extern uint32_t GLEXCOLL_databuf_len;
extern struct TreeLevelInfo _ReduceTreeS[20];
extern int TreeNumber;
extern int TreeLevelN;
extern int _TreeID;
extern int CorePerNuma;
extern int Childn_K;
extern int inter_comm_mode;
extern int INTRA_ALLREDUCE_TYPE;
extern char host_name[MPI_MAX_PROCESSOR_NAME];
extern int *host_ids;

extern int GLEXCOLL_max_PPSDSM_Len;
extern int GLEXCOLL_max_imm_RDMA_Len;
extern int GLEXCOLL_max_RDMA_Len;
extern volatile char *SHM_bufs_Recv[64];
extern volatile char *SHM_bufs_Send[64];
extern volatile int *shared_release_flags[33];
extern volatile int *shared_gather_flags[33];
extern volatile char *broad_cast_buffer;
extern volatile char *Reduce_buffers[33];

int GLEXCOLL_Init(int argc, char **argv);
int GLEXCOLL_Iallreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op);
int GLEXCOLL_Allreduce(void *sendbuf, void *recvbuf, int count,
					   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int GLEXCOLL_AllreduceFinalize();

void one_double_allreduce_intra(void *sendbuf, void *recvbuf, int count);
void one_double_allreduce_intra_REDUCE_to_0(void *sendbuf, void *recvbuf, int count);
void one_double_allreduce_intra_BCAST_from_0(void *sendbuf, void *recvbuf, int count);
int GLEXCOLL_Wait_Iallreduce();
int stringCmp(const void *a, const void *b);
void GLEX_set_affinity(int rank);
void *comm_main_loop(void *pvoid);
int Get_shared_buffer_size();
void GLEXCOLL_Finalize();
void GLEX_ep_init();
void GLEX_ep_destroy();

#define GLEX_MELLOC_MAX 10
extern glex_mem_handle_t *send_mhs[GLEX_MELLOC_MAX];
extern glex_mem_handle_t *recv_mhs[GLEX_MELLOC_MAX];
extern int iphmh_num;
//分配一块用于通信的内存空间。
//该空间分配之后，其它comm中的其它任意进程对该空间可见。
int GLECOLL_Register_SendRecv_Pair(void *sendbuf, uint64_t size1, void *recvbuf, uint64_t size2, MPI_Comm comm);

extern struct _GLEXCOLL_Data _GLEXCOLL;
// #define PJT_NEW_VERSION
#endif
END_C_DECLS