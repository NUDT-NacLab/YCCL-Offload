#include <stdio.h>
#include <stdlib.h>
#include "/usr/local/glex/include/glex.h"
#include "mpi.h"
#include <netdb.h>
#include <fcntl.h> /* For O_* constants */
#include <sys/mman.h>
#include <sys/stat.h> /* For mode constants */
#include <assert.h>
#include <string.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>

extern "C"
{
    int MPI_Init(int *argc, char ***argv);

    int MPI_Bcast(void *data, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

    int MPI_Finalize(void);
}

glex_mem_handle_t local_mh;
glex_ep_handle_t ep_handle_t = NULL;
glex_device_handle_t device_handle_t;

uint64_t *all_ep_addr_value = NULL;
uint64_t ep_addr_value = -1;
uint64_t *all_local_mh_value = NULL;
uint64_t local_mh_value = -1;

uint32_t parent = 0;
uint32_t *child_rank_nomial = NULL;
uint32_t *child_rank_ary = NULL;
uint32_t local_child_cnt_nomial = 0;
uint32_t local_child_cnt_ary = 0;
uint32_t *globle_childNums = NULL;

struct glex_imm_mp_req mp_req;
struct glex_imm_mp_req *bad_mp_req;
struct glex_rdma_req rdma_req;
struct glex_rdma_req *bad_rdma_req;

int flit_size;
int flit_num = 0;
int last_flit = 0;
int offset = 0;
int mess_size = 0;
glex_ep_addr_t rmt_ep_addr;

glex_event_t *event_t;
int buff_size = (1L << 31) - 1; // shared mem
char my_ip[16] = {0};
int init = 0;
int global_rank, global_size; // global rank && size
int rank, size;               // inter rank && size
int node_num = 0;
int type_size = 0;

int nomial_child_cnt = 0;
int ary_child_cnt = 0;

int parent_nomial_init = 0;
int parent_ary_init = 0;

uint32_t len1;

void glex_inter_init(int rank, int size, MPI_Comm comm)
{
    glex_ret_t ret;
    uint32_t num_of_devices;
    ret = glex_num_of_device(&num_of_devices);
    if (ret != GLEX_SUCCESS)
    {
        printf("glex_num_of_device error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }

    ret = glex_open_device(0, &device_handle_t);
    if (ret != GLEX_SUCCESS)
    {
        printf("glex_open_device error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }
    struct glex_ep_attr ep_attr;
    ep_attr.type = GLEX_EP_TYPE_NORMAL;
    ep_attr.mpq_type = GLEX_MPQ_TYPE_NORMAL;
    ep_attr.eq_type = GLEX_EQ_TYPE_NORMAL;
    ep_attr.num = GLEX_ANY_EP_NUM;
    ep_attr.key = 0x23;
    ep_attr.dq_capacity = GLEX_EP_DQ_CAPACITY_DEFAULT;
    ep_attr.mpq_capacity = GLEX_EP_MPQ_CAPACITY_DEFAULT;
    ep_attr.eq_capacity = GLEX_EP_EQ_CAPACITY_DEFAULT;

    ret = glex_create_ep(device_handle_t, &ep_attr, &ep_handle_t);
    if (ret != GLEX_SUCCESS)
    {
        printf("glex_create_ep error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }

    glex_ep_addr_t ep_addr;
    ret = glex_get_ep_addr(ep_handle_t, &ep_addr);
    if (ret != GLEX_SUCCESS)
    {
        printf("glex_get_ep_addr error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }
    ep_addr_value = ep_addr.v;

    all_ep_addr_value = (uint64_t *)calloc(size, sizeof(uint64_t));

    PMPI_Allgather(&ep_addr_value, 1, MPI_UNSIGNED_LONG, all_ep_addr_value, 1,
                   MPI_UNSIGNED_LONG, comm);
}

void glex_inter_mes_init(int bcast_len, void *buffer)
{
    mp_req.data = buffer;
    mp_req.len = bcast_len;
    mp_req.next = NULL;
}

void glex_inter_get_parent_childs_ary(uint32_t rank, uint32_t size, uint32_t group_size, MPI_Comm comm)
{

    child_rank_ary = (uint32_t *)calloc(group_size, sizeof(uint32_t));
    local_child_cnt_ary = 0;

    int p = rank * group_size + 1;
    while (local_child_cnt_ary < group_size && p < size)
    {
        child_rank_ary[local_child_cnt_ary++] = p++;
    }

    globle_childNums = (uint32_t *)calloc(size, sizeof(uint32_t));
    PMPI_Allgather(&local_child_cnt_ary, 1, MPI_UNSIGNED, globle_childNums, 1,
                   MPI_UNSIGNED, comm);
}
void glex_inter_get_parent_childs_nomial(uint32_t rank, uint32_t size, uint32_t group_size, uint32_t ary, MPI_Comm comm)
{
    child_rank_nomial = (uint32_t *)calloc(group_size + ary, sizeof(uint32_t));
    local_child_cnt_nomial = 0;
    if (rank % group_size == 0)
    {
        int pos = rank + 1;
        while (pos < size && local_child_cnt_nomial < group_size - 1)
        {
            child_rank_nomial[local_child_cnt_nomial++] = pos++;
            nomial_child_cnt++;
        }

        int p = ary * rank + group_size;
        int k = 0;
        while (k < ary && p < size)
        {
            child_rank_nomial[local_child_cnt_nomial++] = p;
            p += group_size;
            k++;
            ary_child_cnt++;
        }
        int j = 0;
    }

    globle_childNums = (uint32_t *)calloc(size, sizeof(uint32_t));
    PMPI_Allgather(&local_child_cnt_nomial, 1, MPI_UNSIGNED, globle_childNums, 1,
                   MPI_UNSIGNED, comm);
}

void *tag_shm_addr = NULL;
int intra_size = 0;
void register_mem(int rank, int size, void *buffer, MPI_Comm comm)
{

    glex_ret_t ret;
    // if (intra_size > 1)
    ret = glex_register_mem(ep_handle_t, tag_shm_addr, buff_size,
                            GLEX_MEM_READ | GLEX_MEM_WRITE, &local_mh);
    /* else
        ret = glex_register_mem(ep_handle_t, buffer, 134217728,
                                GLEX_MEM_READ | GLEX_MEM_WRITE, &local_mh); */

    if (ret != GLEX_SUCCESS)
    {
        printf("_register_mem(), return: %d\n", ret);
        PMPI_Finalize();
        exit(1);
    }
    local_mh_value = local_mh.v;
    all_local_mh_value = (uint64_t *)calloc(size, sizeof(uint64_t));
    PMPI_Allgather(&local_mh_value, 1, MPI_UNSIGNED_LONG, all_local_mh_value, 1,
                   MPI_UNSIGNED_LONG, comm);
}
void glex_determine_flit_size(int bcast_len)
{
    if (bcast_len <= 131072)
    {
        flit_size = bcast_len;
    }
    else if (bcast_len <= 8388608)
    {
        flit_size = bcast_len / 8;
    }
    else
    {
        flit_size = bcast_len / 16;
    }
    flit_num = bcast_len / flit_size;

    last_flit = bcast_len % flit_size;

    if ((last_flit) && (last_flit >= flit_size / 2))
        flit_num++;
}
void glex_inter_destroy_ep()
{
    if (ep_handle_t != NULL)
    {
        glex_ret_t ret = glex_destroy_ep(ep_handle_t);
        if (ret != GLEX_SUCCESS)
        {
            printf("glex_destroy_ep error %s.\n", glex_error_str(ret));
            PMPI_Finalize();
        }
    }
}
void glex_inter_deregister_mem()
{
    if (ep_handle_t != NULL)
    {
        glex_ret_t ret = glex_deregister_mem(ep_handle_t, local_mh);
        if (ret != GLEX_SUCCESS)
        {
            printf("glex_deregister_mem error %s.\n", glex_error_str(ret));
            MPI_Finalize();
        }
    }
}
void inter_bcast_rdma(MPI_Comm comm, int len, void *bcast_buffer)
{

    // int rank, size;
    PMPI_Comm_rank(comm, &rank);
    PMPI_Comm_size(comm, &size);
    if (init == 0)
    {
        glex_inter_init(rank, size, comm);

        init = 1;
    }

    if (parent_nomial_init == 0)
    {
        register_mem(rank, size, bcast_buffer, comm);
        glex_inter_get_parent_childs_nomial(rank, size, 3, 2, comm);
        parent_nomial_init = 1;
    }
    glex_determine_flit_size(len);
    glex_ret_t ret;

    if (rank == 0)
    {
        // if (intra_size > 1)
        memcpy((void *)tag_shm_addr, bcast_buffer, len);
        int j = 0;
        int a = 0;
        rdma_req.local_mh.v = all_local_mh_value[rank];
        rdma_req.type = GLEX_RDMA_TYPE_PUT;
        rdma_req.rmt_key = 0x23;
        rdma_req.next = NULL;
        rdma_req.rmt_evt.cookie_0 = 0;
        rdma_req.rmt_evt.cookie_1 = 4;
        for (a = 0; a < flit_num; ++a)
        {

            if ((a == (flit_num - 1)) && (last_flit != 0) && (last_flit >= flit_size / 2))
            {
                rdma_req.len = last_flit;
                rdma_req.local_offset = a * flit_size;
                rdma_req.rmt_offset = a * flit_size;
            }
            else if ((a == (flit_num - 1)) && (last_flit != 0) && (last_flit < flit_size / 2))
            {
                rdma_req.len = flit_size + last_flit;
                rdma_req.local_offset = a * flit_size;
                rdma_req.rmt_offset = a * flit_size;
            }
            else
            {
                rdma_req.len = flit_size;
                rdma_req.local_offset = a * flit_size;
                rdma_req.rmt_offset = a * flit_size;
            }
            for (j = nomial_child_cnt; j < local_child_cnt_nomial; ++j)
            {
                rdma_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_nomial[j]];
                rdma_req.rmt_mh.v = all_local_mh_value[child_rank_nomial[j]];
                rdma_req.flag = GLEX_FLAG_COLL;
                rdma_req.flag |= GLEX_FLAG_NO_CONNECTION;
                if (a == flit_num - 1)
                {
                    rdma_req.flag |= GLEX_FLAG_REMOTE_EVT;
                }
                if (globle_childNums[child_rank_nomial[j]])
                {
                    rdma_req.flag |= GLEX_FLAG_COLL_COUNTER;
                }
                if (j == nomial_child_cnt)
                {
                    rdma_req.coll_counter = 0;
                    rdma_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                }
                if (j == local_child_cnt_nomial - 1)
                {
                    rdma_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;
                }

                while ((ret = glex_rdma(ep_handle_t, &rdma_req, &bad_rdma_req)) ==
                       GLEX_BUSY)
                    ;
            }
        }

        for (j = 0; j < nomial_child_cnt; ++j)
        {

            rdma_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_nomial[j]];
            rdma_req.rmt_mh.v = all_local_mh_value[child_rank_nomial[j]];

            rdma_req.flag = GLEX_FLAG_COLL | GLEX_FLAG_REMOTE_EVT;
            rdma_req.flag |= GLEX_FLAG_NO_CONNECTION;
            if (globle_childNums[child_rank_nomial[j]])
            {
                rdma_req.flag |= GLEX_FLAG_COLL_COUNTER;
            }
            rdma_req.len = len;
            rdma_req.local_offset = 0;
            rdma_req.rmt_offset = 0;
            if (j == 0)
            {
                rdma_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                rdma_req.coll_counter = 0;
            }
            if (j == ary_child_cnt - 1)
            {
                rdma_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;
            }

            while ((ret = glex_rdma(ep_handle_t, &rdma_req, &bad_rdma_req)) ==
                   GLEX_BUSY)
                ;
        }
    }
    else
    {

        if (local_child_cnt_nomial)
        {

            rdma_req.local_mh.v = all_local_mh_value[rank];
            rdma_req.type = GLEX_RDMA_TYPE_PUT;
            rdma_req.rmt_key = 0x23;
            rdma_req.next = NULL;
            rdma_req.rmt_evt.cookie_0 = 0;
            rdma_req.rmt_evt.cookie_1 = 4;
            int j = 0;
            int a = 0;
            if (ary_child_cnt)
            {
                for (a = 0; a < flit_num; ++a)
                {
                    if ((a == (flit_num - 1)) && (last_flit != 0) && (last_flit >= flit_size / 2))
                    {
                        rdma_req.len = last_flit;
                        rdma_req.local_offset = a * flit_size;
                        rdma_req.rmt_offset = a * flit_size;
                    }
                    else if ((a == (flit_num - 1)) && (last_flit != 0) && (last_flit < flit_size / 2))
                    {
                        rdma_req.len = flit_size + last_flit;
                        rdma_req.local_offset = a * flit_size;
                        rdma_req.rmt_offset = a * flit_size;
                    }
                    else
                    {
                        rdma_req.len = flit_size;
                        rdma_req.local_offset = a * flit_size;
                        rdma_req.rmt_offset = a * flit_size;
                    }
                    for (j = nomial_child_cnt; j < local_child_cnt_nomial; ++j)
                    {

                        rdma_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_nomial[j]];
                        rdma_req.rmt_mh.v = all_local_mh_value[child_rank_nomial[j]];

                        rdma_req.flag = GLEX_FLAG_COLL;
                        // rdma_req.flag |= GLEX_FLAG_NO_CONNECTION;
                        if (a == flit_num - 1)
                        {
                            rdma_req.flag |= GLEX_FLAG_REMOTE_EVT;
                        }
                        if (globle_childNums[child_rank_nomial[j]])
                        {
                            rdma_req.flag |= GLEX_FLAG_COLL_COUNTER;
                        }
                        if (j == nomial_child_cnt)
                        {
                            rdma_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                            rdma_req.coll_counter = 1;
                        }
                        if (j == local_child_cnt_nomial - 1)
                        {
                            rdma_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;
                        }

                        while ((ret = glex_rdma(ep_handle_t, &rdma_req, &bad_rdma_req)) ==
                               GLEX_BUSY)
                            ;
                    }
                }

                for (j = 0; j < nomial_child_cnt; ++j)
                {

                    rdma_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_nomial[j]];
                    rdma_req.rmt_mh.v = all_local_mh_value[child_rank_nomial[j]];
                    rdma_req.flag = GLEX_FLAG_COLL | GLEX_FLAG_REMOTE_EVT | GLEX_FLAG_FENCE;
                    // rdma_req.flag |= GLEX_FLAG_NO_CONNECTION;
                    rdma_req.len = len;
                    rdma_req.local_offset = 0;
                    rdma_req.rmt_offset = 0;
                    if (j == 0)
                    {
                        rdma_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                        rdma_req.coll_counter = 0;
                    }
                    if (j == ary_child_cnt - 1)
                    {
                        rdma_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;
                    }

                    while ((ret = glex_rdma(ep_handle_t, &rdma_req, &bad_rdma_req)) ==
                           GLEX_BUSY)
                        ;
                }
            }
            else
            {
                for (a = 0; a < flit_num; ++a)
                {
                    if ((a == (flit_num - 1)) && (last_flit != 0) && (last_flit >= flit_size / 2))
                    {
                        rdma_req.len = last_flit;
                        rdma_req.local_offset = a * flit_size;
                        rdma_req.rmt_offset = a * flit_size;
                    }
                    else if ((a == (flit_num - 1)) && (last_flit != 0) && (last_flit < flit_size / 2))
                    {
                        rdma_req.len = flit_size + last_flit;
                        rdma_req.local_offset = a * flit_size;
                        rdma_req.rmt_offset = a * flit_size;
                    }
                    else
                    {
                        rdma_req.len = flit_size;
                        rdma_req.local_offset = a * flit_size;
                        rdma_req.rmt_offset = a * flit_size;
                    }

                    for (j = 0; j < nomial_child_cnt; ++j)
                    {

                        rdma_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_nomial[j]];
                        rdma_req.rmt_mh.v = all_local_mh_value[child_rank_nomial[j]];

                        rdma_req.flag = GLEX_FLAG_COLL;
                        // rdma_req.flag |= GLEX_FLAG_NO_CONNECTION;
                        if (a == flit_num - 1)
                        {
                            rdma_req.flag |= GLEX_FLAG_REMOTE_EVT;
                        }
                        if (globle_childNums[child_rank_nomial[j]])
                        {
                            rdma_req.flag |= GLEX_FLAG_COLL_COUNTER;
                        }

                        if (j == 0)
                        {
                            rdma_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                            rdma_req.coll_counter = 1;
                        }
                        if (j == ary_child_cnt - 1)
                        {
                            rdma_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;
                        }

                        while ((ret = glex_rdma(ep_handle_t, &rdma_req, &bad_rdma_req)) ==
                               GLEX_BUSY)
                            ;
                    }
                }
            }

            while ((ret = glex_probe_next_event(ep_handle_t, &event_t)) ==
                   GLEX_NO_EVENT)
                ;
            glex_discard_probed_event(ep_handle_t);
        }
        else
        {
            while ((ret = glex_probe_next_event(ep_handle_t, &event_t)) ==
                   GLEX_NO_EVENT)
                ;
            glex_discard_probed_event(ep_handle_t);
        }
    }
    // glex_inter_deregister_mem();
}

void inter_bcast_mp(MPI_Comm comm, int len, void *bcast_buffer)
{
    PMPI_Comm_rank(comm, &rank);
    PMPI_Comm_size(comm, &size);

    if (init == 0)
    {
        glex_inter_init(rank, size, comm);
        init = 1;
    }
    if (parent_ary_init == 0)
    {
        glex_inter_get_parent_childs_ary(rank, size, 16, comm);
        parent_ary_init = 1;
    }

    glex_inter_mes_init(len, bcast_buffer);

    glex_ret_t ret;
    if (rank == 0)
    {
        int j = 0;
        for (j = 0; j < local_child_cnt_ary; j++)
        {
            mp_req.flag = GLEX_FLAG_COLL | GLEX_FLAG_COLL_CP_DATA;
            mp_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_ary[j]];
            mp_req.flag |= GLEX_FLAG_NO_CONNECTION;
            /* if ((root_inter_node != 0) && (child_rank_ary[j] == root_inter_node))
            {
                mp_req.rmt_ep_addr.v = all_ep_addr_value[0];
            }
            else
            {
                mp_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_ary[j]];
            } */
            if (j == 0)
            {
                mp_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                mp_req.coll_counter = 0;
            }

            if (j == (local_child_cnt_ary - 1))
            {
                mp_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;
            }

            if (globle_childNums[child_rank_ary[j]])
            {
                mp_req.flag |= GLEX_FLAG_COLL_CP_SWAP_DATA |
                               GLEX_FLAG_COLL_COUNTER;
            }

            while ((ret = glex_send_imm_mp(ep_handle_t, &mp_req, &bad_mp_req)) ==
                   GLEX_BUSY)
                ;
        }
    }
    else
    {

        if (local_child_cnt_ary)
        {
            int j = 0;
            for (j = 0; j < local_child_cnt_ary; j++)
            {
                mp_req.flag = GLEX_FLAG_COLL | GLEX_FLAG_COLL_CP_DATA | GLEX_FLAG_COLL_SWAP;
                mp_req.flag |= GLEX_FLAG_NO_CONNECTION;
                mp_req.rmt_ep_addr.v = all_ep_addr_value[child_rank_ary[j]];

                if (j == 0)
                {
                    mp_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                    mp_req.coll_counter = 1;
                }

                if (j == (local_child_cnt_ary - 1))
                {
                    mp_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;
                }

                if (globle_childNums[child_rank_ary[j]])
                {
                    mp_req.flag |= GLEX_FLAG_COLL_CP_SWAP_DATA |
                                   GLEX_FLAG_COLL_COUNTER;
                }

                while ((ret = glex_send_imm_mp(ep_handle_t, &mp_req, &bad_mp_req)) == GLEX_BUSY)
                    ;
            }
            if (intra_size > 1)
            {
                while ((ret = glex_receive_mp(ep_handle_t, 0, &rmt_ep_addr, (void *)tag_shm_addr,
                                              &len1)) == GLEX_NO_MP)
                {
                }
            }
            else
            {
                while ((ret = glex_receive_mp(ep_handle_t, 0, &rmt_ep_addr, bcast_buffer,
                                              &len1)) == GLEX_NO_MP)
                {
                }
            }
        }
        else
        {
            if (intra_size > 1)
            {
                while ((ret = glex_receive_mp(ep_handle_t, 0, &rmt_ep_addr, (void *)tag_shm_addr,
                                              &len1)) == GLEX_NO_MP)
                {
                }
            }
            else
            {
                while ((ret = glex_receive_mp(ep_handle_t, 0, &rmt_ep_addr, bcast_buffer,
                                              &len1)) == GLEX_NO_MP)
                {
                }
            }
        }
    }
}

int *leader_rank = NULL;
int leader_count = 0;

int tag_shm_fd = 0;

int yhccl_new_rank = 0;

void get_leader_and_shm_init(int rank, int size, int count, MPI_Datatype datatype)
{

    char my_ip[16] = {0};
    char str_host_name[100];
    gethostname(str_host_name, 100);
    inet_ntop(AF_INET, gethostbyname(str_host_name)->h_addr_list[0], my_ip, 16);

    char *all_ip = (char *)calloc(size, 16);
    PMPI_Allgather(my_ip, 16, MPI_CHAR, all_ip, 16, MPI_CHAR, MPI_COMM_WORLD);

    leader_rank = (int *)calloc(size, sizeof(int));
    leader_count = 0;
    yhccl_new_rank = 0;
    intra_size = 0;
    node_num = 0;

    char *tp = all_ip;
    int i = 0;
    for (i = 0; i < size; ++i, tp += 16)
    {
        if (strcmp(tp, my_ip) == 0)
        {
            ++intra_size;
            if (i < rank)
                ++yhccl_new_rank;
        }

        if ((i <= rank) && (i != 0) && (strcmp(tp, tp - 16) != 0))
        {
            node_num++;
        }

        if (i != 0)
        {
            if (strcmp(tp, tp - 16) == 0)
            {

                continue;
            }
        }

        leader_rank[leader_count++] = i;
    }

    tag_shm_fd = shm_open(my_ip, O_CREAT | O_RDWR, 0777);
    assert(tag_shm_fd >= 0);

    if (ftruncate(tag_shm_fd, buff_size) == -1)
    {
        close(tag_shm_fd);
        tag_shm_fd = 0;
        assert(0);
    }

    tag_shm_addr = (char *)mmap(NULL, buff_size, PROT_READ | PROT_WRITE, MAP_SHARED, tag_shm_fd, 0);
    if (tag_shm_addr == MAP_FAILED)
    {
        close(tag_shm_fd);
        tag_shm_fd = 0;
        assert(0);
    }
}

void destory_shm_fd()
{
    munmap(tag_shm_addr, buff_size);
    tag_shm_addr = 0;

    if (tag_shm_fd != 0)
        close(tag_shm_fd);
    tag_shm_fd = 0;

    char my_ip[16] = {0};
    char str_host_name[100];
    gethostname(str_host_name, 100);

    inet_ntop(AF_INET, gethostbyname(str_host_name)->h_addr_list[0], my_ip, 16);
    shm_unlink(my_ip);
}

int is_init_intra = 0;
MPI_Comm yhccl_inter_comm;
MPI_Comm yhccl_intra_comm;
int yhccl_inter_comm_split_flag = 0;
int yhccl_intra_comm_split_flag = 0;

void intra_bcast_rdma(void *buffer, int bcast_len)
{

    /* if (global_rank == 0)
    {
        memcpy((void *)tag_shm_addr, buffer, bcast_len);
    } */
    PMPI_Barrier(yhccl_intra_comm);
    if (global_rank > 0)
    {
        memcpy(buffer, (void *)tag_shm_addr, bcast_len);
    }
}

void intra_bcast_mp(void *buffer, int bcast_len)
{
    if (global_rank == 0)
    {
        memcpy((void *)tag_shm_addr, buffer, bcast_len);
    }
    PMPI_Barrier(yhccl_intra_comm);
    if (global_rank > 0)
    {
        memcpy(buffer, tag_shm_addr, bcast_len);
    }
}
/********************************intra*****************************/
void bcast_rdma(void *buffer, int count, int root, MPI_Comm comm, MPI_Datatype datatype)
{

    if (yhccl_new_rank == 0)
    {
        inter_bcast_rdma(yhccl_inter_comm, mess_size, buffer);
    }
    // if (intra_size > 1)
    intra_bcast_rdma(buffer, mess_size);
}

void bcast_mp(void *buffer, int count, int root, MPI_Comm comm, MPI_Datatype datatype)
{
    if (yhccl_new_rank == 0)
    {
        inter_bcast_mp(yhccl_inter_comm, mess_size, buffer);
    }
    if (intra_size > 1)
        intra_bcast_mp(buffer, mess_size);
}

int MPI_Init(int *argc, char ***argv)
{
    return PMPI_Init(argc, argv);
}

int MPI_Bcast(void *data, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
    PMPI_Comm_rank(comm, &global_rank);
    PMPI_Comm_size(comm, &global_size);
    PMPI_Type_size(datatype, &type_size);
    mess_size = count * type_size;
    // printf("mess_size is %d \n", mess_size);
    if (root != 0)
    {
        return PMPI_Bcast(data, count, datatype, root, comm);
    }

    if (is_init_intra == 0)
    {
        get_leader_and_shm_init(global_rank, global_size, count, datatype);
        is_init_intra = 1;
    }

    if (yhccl_inter_comm_split_flag == 0)
    {
        int *color = (int *)calloc(global_size, sizeof(int));
        int i = 0;
        for (; i < leader_count; i++)
        {
            color[leader_rank[i]] = 1;
        }
        PMPI_Comm_split(comm, color[global_rank], global_rank, &yhccl_inter_comm);
        yhccl_inter_comm_split_flag = 1;
    }
    // if ((intra_size > 1) && (yhccl_intra_comm_split_flag == 0))

    if (yhccl_intra_comm_split_flag == 0)
    {

        PMPI_Comm_split(comm, node_num, global_rank, &yhccl_intra_comm);
        yhccl_intra_comm_split_flag = 1;
    }
    if (mess_size <= 112) // RDMA
    {
        // printf("here1\n");
        bcast_mp(data, count, root, comm, datatype);
    }
    else
    {
        bcast_rdma(data, count, root, comm, datatype);
    }
    return 0;
}
int MPI_Finalize(void)
{
    if (mess_size > 112)
        glex_inter_deregister_mem();
    destory_shm_fd();
    glex_inter_destroy_ep();
    return PMPI_Finalize();
}
