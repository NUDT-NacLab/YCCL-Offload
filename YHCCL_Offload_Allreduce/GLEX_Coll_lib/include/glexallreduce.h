#ifndef GLEXALLREDUCE_H
#define GLEXALLREDUCE_H

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

enum SoftWare_Allreduce_Algorithm_type
{
    Recursize_doubling_Slicing,
    Recursize_doubling_OMP,
    K_nomial_tree_OMP,
    Self_tuned_two_dimentional,
    Small_message_node_aware_POSIX_RDMA
};

enum Bcast_Algorithm
{
    K_ary_broadcast,
    K_ary_broadcast_pipeline
};
extern void GLEXCOLL_Allreduce_NonOffload(void *sendbuf, int count, void *recvbuf);
extern void GLEXCOLL_Bcast(
    void *data_p,
    int count,
    MPI_Datatype datatype,
    int source_proc,
    MPI_Comm comm);
void GLEXCOLL_InitAllreduce();
extern char *allreduce_sendbuf, *allreduce_recvbuf;
extern int allreduce_recursive_doubling_stepn;
extern int allreduce_send_recv_pair;

extern void GLEX_Small_message_bcast_double(double *sendbuf, double *recvbuf, int count);
extern void GLEX_Small_message_reduce_double_sum(double *sendbuf, double *recvbuf, int count);
extern void medium_reduce(double *sendbuf, double *recvbuf, int count);
extern void medium_bcast(double *sendbuf, double *recvbuf, int count);

extern unsigned long allreduce_buffer_size_total;
extern unsigned long allreduce_buffer_size_single;

struct sync_lock
{
    volatile int lock0[16];
};
#define allreduce_intra_node_locks_n 64
typedef struct sync_header
{
    volatile struct sync_lock mlocks[allreduce_intra_node_locks_n];
} allreduce_Header;

extern int SoftWare_Allreduce_Algorithm;
extern volatile allreduce_Header *allreduce_shm_flags;
extern void *allreduce_shm_buffer_starts[];
END_C_DECLS

#endif