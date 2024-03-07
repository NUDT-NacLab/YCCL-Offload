#ifndef GLEX_ALLTOALL_H
#define GLEX_ALLTOALL_H
#include <stdio.h>
#include <mpi.h>


#ifdef __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else /* !__cplusplus */
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif /* __cplusplus */

BEGIN_C_DECLS

//alltoall算法选择变量，通过环境变量 ALLTOALL_TYPE 导入
enum ALGORITHM {BRUCK,BRUCK_RDMA,
                DIRECT,
                DIRECT_NODE_AWARE,
                DIRECT_Kleader_NODE_AWARE,
                DIRECT_Kleader_NODE_AWARE_RDMA,
                DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE,
                XOR_EXCHANGE_RDMA_ONGOING,
                MPI_LINEAR_EXCHANGE_COMPRESS,
                TPDS17_Cache_oblivious_intra_node,
                shared_memory_direct_intra_node,
                TPDS17_Cache_oblivious_intra_node_NUMA,
                L_a2a,
                ONMPML,
                MPML,
                NMPML,
                SONMPML};
extern int Alltoall_algorithm;

void glexcoll_InitAlltoall();
void glexcoll_init_alltoall_shared_memory_buffer();
extern void *get_is_senddata_buffer(int i);
extern void *get_is_recvdata_buffer(int i);
extern void *tpds17_get_is_senddata_buffer(int i);
extern void *tpds17_get_is_recvdata_buffer(int i);

int GLEXCOLL_Alltoall(void *sendbuf,
                      int sendsize,
                      void *recvbuf,
                      int recvsize,
                      MPI_Comm comm,
                      MPI_Datatype type);
                
void GLEXCOLL_AlltoallFinalize();

//运行非MPI_Comm_world的alltoallv
//目前支持没节点单个进程
void GLEXCOLL_InitAlltoallV(MPI_Comm comm);
void  GLEXCOLL_AlltoallVFinalize();
void GLEXCOLL_Alltoallv(void *sendbuf,
                        uint64_t sendcounts[],
                        uint64_t sdispls[],
                        void *recvbuf,
                        uint64_t recvcounts[],
                        uint64_t rdispls[],
                        MPI_Datatype  type);

void GLEXCOLL_Alltoall_new(void *sendbuf,
                      int sendsize,
                      void *recvbuf,
                      int recvsize,
                      MPI_Comm comm,
                      MPI_Datatype type);
void glexcoll_InitAlltoall_new(MPI_Comm comm);
void GLEXCOLL_Alltoall_Finalize_new();

enum AlltoallvAlgorithmTYPE
{
    LINEAR_SHIFT,
    XOR,
    OPT_BANDWIDTH,
    LINEAR_SHIFT_REORDER,
    RDMA_LINEAR_SHIFT,
    RDMA_LINEAR_SHIFT_REORDER,
    LINEAR_SHIFT_COMPRESS
};
struct GLEXCOLL_a2a_bufmh{
    uint64_t sendvec[64];
    uint64_t recvvec[64];
    glex_mem_handle_t sendmh;
    glex_mem_handle_t *alltoall_mem_handle_vec;
    glex_ep_addr_t *alltoall_ep_addrs;
    void * sendbuf;
    void * recvbuf;
};
void glexcoll_register_alltoall_buffer(void *sendbuf, void *recvbuf,struct GLEXCOLL_a2a_bufmh *bufmh);
void glexcoll_register_alltoall_buffer_new(void * sendbuf,void *recvbuf,int size, MPI_Comm comm,struct GLEXCOLL_a2a_bufmh *bufmh);

void GLEXCOLL_Alltoall_pjt(struct GLEXCOLL_a2a_bufmh *bufmh, int size);
void GLEXCOLL_Alltoallv_Performance_detection(void *sendbuf,
                                                          uint64_t sendcounts[],
                                                          uint64_t sdispls[],
                                                          void *recvbuf,
                                                          uint64_t recvcounts[],
                                                          uint64_t rdispls[]);
void GLEXCOLL_Alltoallv_RDMA_REORDER(int mh_pair, void *sendbuf,void *recvbuf,uint64_t sendcounts[],uint64_t sdispls[]);
extern void write_out_print_alltoallv_performance_data(int loopn,double *Step_max_sum);
extern void GLEXCOLL_Alltollv_RDMA(int mh_pair, void *sendbuf,void *recvbuf,uint64_t sendcounts[],uint64_t sdispls[]);
void GLEXCOLL_Alltoallv_regist_recv_disps(uint64_t * rdispls);
extern int alltoallv_single_data_size;
END_C_DECLS
#endif