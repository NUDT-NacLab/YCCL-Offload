#include <stdio.h>
#include <stdlib.h>
#include "/usr/local/glex/include/glex.h"
#include "mpi.h"
// #include "/usr/local/ompi-x/include/mpi.h"

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
//  ↓↓↓↓↓ ######### inter ######### ↓↓↓↓↓
glex_ep_handle_t ep_handle_t = NULL;

uint64_t* all_ep_addr_value = NULL;
uint64_t ep_addr_value = -1;  // ep值

uint32_t parent = 0;                // 父亲
uint32_t* child_rank = NULL;        // 儿子们
uint32_t local_child_cnt = 0;       // 本进程儿子们的数量
uint32_t* globle_childNums = NULL;  // 全局每个进程的孩子数量

struct glex_imm_mp_req mp_req;
struct glex_imm_mp_req* bad_mp_req;
uint64_t data[112 / 8];
uint32_t len;
glex_ep_addr_t rmt_ep_addr;
//  ↑↑↑↑↑ ######### inter ######### ↑↑↑↑↑

extern "C" {
int MPI_Init(int* argc, char*** argv);

int MPI_Barrier(MPI_Comm comm);

int MPI_Finalize(void);
}

int init = 0;

void glex_inter_get_parent_childs(uint32_t rank, uint32_t size, uint32_t group_size, MPI_Comm comm) {
    parent = (rank - 1) / group_size;                              // 0号无需管
    child_rank = (uint32_t*)calloc(group_size, sizeof(uint32_t));  // 孩子最多叉度的数量
    local_child_cnt = 0;

    int p = rank * group_size + 1;
    while (local_child_cnt < group_size && p < size) {
        child_rank[local_child_cnt++] = p++;
    }

    // 获取每个leader的儿子数量
    globle_childNums = (uint32_t*)calloc(size, sizeof(uint32_t));
    PMPI_Allgather(&local_child_cnt, 1, MPI_UNSIGNED, globle_childNums, 1, MPI_UNSIGNED, comm);
}

// 每个节点的leader执行
void glex_inter_init(int rank, int size, MPI_Comm comm) {
    glex_ret_t ret;
    uint32_t num_of_devices;
    ret = glex_num_of_device(&num_of_devices);
    if (ret != GLEX_SUCCESS) {
        printf("glex_num_of_device error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }

    // 获取设备device_handle
    glex_device_handle_t device_handle_t;
    ret = glex_open_device(0, &device_handle_t);
    if (ret != GLEX_SUCCESS) {
        printf("glex_open_device error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }

    // 准备端点信息
    struct glex_ep_attr ep_attr;
    // ep_attr.type = GLEX_EP_TYPE_FAST;
    ep_attr.type = GLEX_EP_TYPE_NORMAL;
    ep_attr.mpq_type = GLEX_MPQ_TYPE_NORMAL;
    ep_attr.eq_type = GLEX_EQ_TYPE_NORMAL;
    ep_attr.num = GLEX_ANY_EP_NUM;
    ep_attr.key = GLEX_EP_MAKE_KEY(999);
    ep_attr.dq_capacity = GLEX_EP_DQ_CAPACITY_DEFAULT;
    ep_attr.mpq_capacity = GLEX_EP_MPQ_CAPACITY_DEFAULT;
    ep_attr.eq_capacity = GLEX_EP_EQ_CAPACITY_DEFAULT;
    // 创建端点
    ret = glex_create_ep(device_handle_t, &ep_attr, &ep_handle_t);
    if (ret != GLEX_SUCCESS) {
        printf("glex_create_ep error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }

    // 获取端点信息
    glex_ep_addr_t ep_addr;
    ret = glex_get_ep_addr(ep_handle_t, &ep_addr);
    if (ret != GLEX_SUCCESS) {
        printf("glex_get_ep_addr error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }
    ep_addr_value = ep_addr.v;

    all_ep_addr_value = (uint64_t*)calloc(size, sizeof(uint64_t));
    // 全局leader规约  -- ep_addr_value
    PMPI_Allgather(&ep_addr_value, 1, MPI_UNSIGNED_LONG, all_ep_addr_value, 1, MPI_UNSIGNED_LONG, comm);
}



// leader执行
void glex_inter_destroy_ep() {
    if (ep_handle_t != NULL) {
        glex_ret_t ret = glex_destroy_ep(ep_handle_t);
        if (ret != GLEX_SUCCESS) {
            printf("glex_destroy_ep error %s.\n", glex_error_str(ret));
            PMPI_Finalize();
        }
    }
}

void glex_inter_mes_init() {
    mp_req.data = NULL;
    mp_req.len = 0;
    mp_req.next = NULL;
}

void glex_inter_barrier(MPI_Comm comm) {
    int rank, size;
    PMPI_Comm_rank(comm, &rank);
    PMPI_Comm_size(comm, &size);

    if (init == 0) {
        glex_inter_init(rank, size, comm);
        // 获取leader 拓扑父子信息
        glex_inter_get_parent_childs(rank, size, 8, comm);
        glex_inter_mes_init();
        init = 1;
    }

    glex_ret_t ret;
    if (rank == 0) {
        int i = 0;
        for (i = local_child_cnt - 1; i >= -1; --i) {
            mp_req.flag = GLEX_FLAG_COLL;
            mp_req.coll_counter = 0;

            if (i == local_child_cnt - 1) {
                mp_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                mp_req.coll_counter = local_child_cnt;
            }

            // 自己
            if (i == -1) {
                mp_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;  // 最后一个发送到mpq
                mp_req.rmt_ep_addr.v = all_ep_addr_value[rank];
            } else {
                mp_req.rmt_ep_addr.v = all_ep_addr_value[child_rank[i]];
            }

            if (i == -1 || !globle_childNums[child_rank[i]]) {
                // 对方无子结点 或者是本身
                mp_req.flag |= GLEX_FLAG_COLL_CP_DATA;
            } else {
                // 需要触发对方的counter
                mp_req.flag |= GLEX_FLAG_COLL_COUNTER;
            }

            mp_req.flag |= GLEX_FLAG_NO_CONNECTION;
            while ((ret = glex_send_imm_mp(ep_handle_t, &mp_req,
                                           &bad_mp_req)) == GLEX_BUSY)
                ;
        }

    } else {
        if (local_child_cnt) {
            mp_req.rmt_ep_addr.v = all_ep_addr_value[parent];
            mp_req.flag = GLEX_FLAG_COLL_SEQ_START | GLEX_FLAG_COLL_SEQ_TAIL |
                          GLEX_FLAG_COLL | GLEX_FLAG_COLL_COUNTER;
            mp_req.coll_counter = local_child_cnt;
            mp_req.flag |= GLEX_FLAG_NO_CONNECTION;
            while ((ret = glex_send_imm_mp(ep_handle_t, &mp_req,
                                           &bad_mp_req)) == GLEX_BUSY)
                ;
            int i = 0;
            for (i = local_child_cnt - 1; i >= -1; --i) {
                mp_req.flag = GLEX_FLAG_COLL;
                mp_req.coll_counter = 0;

                if (i == local_child_cnt - 1) {
                    mp_req.flag |= GLEX_FLAG_COLL_SEQ_START;
                    mp_req.coll_counter = 1;
                }

                // 自己
                if (i == -1) {
                    mp_req.flag |= GLEX_FLAG_COLL_SEQ_TAIL;  // 最后一个发送到mpq
                    mp_req.rmt_ep_addr.v = all_ep_addr_value[rank];
                } else {
                    mp_req.rmt_ep_addr.v = all_ep_addr_value[child_rank[i]];
                }

                if (i == -1 || !globle_childNums[child_rank[i]]) {
                    // 对方无子结点 或者是本身
                    mp_req.flag |= GLEX_FLAG_COLL_CP_DATA;
                } else {
                    // 需要触发对方的counter
                    mp_req.flag |= GLEX_FLAG_COLL_COUNTER;
                }
                mp_req.flag |= GLEX_FLAG_NO_CONNECTION;
                while ((ret = glex_send_imm_mp(ep_handle_t, &mp_req,
                                               &bad_mp_req)) == GLEX_BUSY)
                    ;
                //         printf("leader%u: child:%u\n", rank,
                //         getRank(all_epaddr, mp_req.rmt_ep_addr.v));
            }

        } else {
            // 叶子节点
            mp_req.rmt_ep_addr.v = all_ep_addr_value[parent];
            mp_req.flag = GLEX_FLAG_COLL | GLEX_FLAG_COLL_COUNTER;
            mp_req.flag |= GLEX_FLAG_NO_CONNECTION;
            while ((ret = glex_send_imm_mp(ep_handle_t, &mp_req,
                                           &bad_mp_req)) == GLEX_BUSY)
                ;
        }
    }

    // while ((ret = glex_receive_mp(ep_handle_t, 0, &rmt_ep_addr, data, &len)) ==
    //        GLEX_NO_MP) {
    // }

    while ((ret = glex_probe_next_mp(ep_handle_t, &rmt_ep_addr, (void**)&data, &len)) ==
           GLEX_NO_MP) {
    }
    static int count_236;
    count_236++;

    if (count_236 % 128 == 0) {
        count_236 = 0;
        if (glex_discard_probed_mp(ep_handle_t) != GLEX_SUCCESS) {
            puts("error glex_discard_probed_mp");
        }
    }

    if (ret != GLEX_SUCCESS) {
        printf("glex_receive_mp error %s.\n", glex_error_str(ret));
        PMPI_Finalize();
    }
}

int* leader_rank = NULL;
int leader_count = 0;

int yhccl_shm_fd = 0;
char* yhccl_shm_addr = NULL;
int yhccl_buf_size = 0;

int yhccl_new_rank = 0;  // 节点内的rank
int yhccl_new_size = 0;  // 节点内的size

// 接收全局comm的size
void get_leader_and_shm_init(int rank, int size) {
    char my_ip[16] = {0};
    char str_host_name[100];
    gethostname(str_host_name, 100);
    inet_ntop(AF_INET, gethostbyname(str_host_name)->h_addr_list[0], my_ip, 16);

    char* all_ip = (char*)calloc(size, 16);
    PMPI_Allgather(my_ip, 16, MPI_CHAR, all_ip, 16, MPI_CHAR, MPI_COMM_WORLD);

    leader_rank = (int*)calloc(size, sizeof(int));
    leader_count = 0;
    yhccl_new_rank = 0;
    yhccl_new_size = 0;

    char* tp = all_ip;
    int i = 0;
    for (i = 0; i < size; ++i, tp += 16) {
        if (strcmp(tp, my_ip) == 0) {
            ++yhccl_new_size;
            if (i < rank)
                ++yhccl_new_rank;
        }

        if (i != 0) {
            if (strcmp(tp, tp - 16) == 0)
                continue;
        }
        leader_rank[leader_count++] = i;
    }

    // 创建共享文件
    yhccl_shm_fd = shm_open(my_ip, O_CREAT | O_RDWR, 0777);
    assert(yhccl_shm_fd >= 0);

    // global
    yhccl_buf_size = sizeof(char) * yhccl_new_size * 64;

    // 申请大小
    if (ftruncate(yhccl_shm_fd, yhccl_buf_size) == -1) {
        close(yhccl_shm_fd);
        yhccl_shm_fd = 0;
        assert(0);
    }

    // 进行映射
    yhccl_shm_addr = (char*)mmap(NULL, yhccl_buf_size, PROT_READ | PROT_WRITE,
                                 MAP_SHARED, yhccl_shm_fd, 0);
    if (yhccl_shm_addr == MAP_FAILED) {
        close(yhccl_shm_fd);
        yhccl_shm_fd = 0;
        assert(0);
    }

    PMPI_Barrier(MPI_COMM_WORLD);
}

void destory_shm_fd() {
    munmap(yhccl_shm_addr, yhccl_buf_size);
    yhccl_buf_size = 0;
    yhccl_shm_addr = 0;

    if (yhccl_shm_fd != 0)
        close(yhccl_shm_fd);
    yhccl_shm_fd = 0;

    char my_ip[16] = {0};
    char str_host_name[100];
    gethostname(str_host_name, 100);

    inet_ntop(AF_INET, gethostbyname(str_host_name)->h_addr_list[0], my_ip, 16);
    shm_unlink(my_ip);
}

// void rank_0_detect_all_1() __attribute__((optimize("O0")));
// void rank_0_detect_all_1() {
//     if (0 == yhccl_new_rank) {
//         int i = 0;
//         for (i = 0; i < yhccl_new_size; ++i) {
//             if (0 == yhccl_shm_addr[i]) {
//                 i -= 1;
//             }
//         }
//     }
// }

// void my_wait() __attribute__((optimize("O0")));
// void my_wait() {
//     //节点间完毕,newrank0出来，告知节点内的可以离开了
//     if (yhccl_new_rank == 0) {
//         memset(yhccl_shm_addr, 0, sizeof(char) * yhccl_new_size);
//     }

//     if (yhccl_new_size > 1) {
//         while (yhccl_shm_addr[yhccl_new_rank] != 0) {
//         }
//     }
// }

// [[gnu::optimize(0)]]
void while_wait() __attribute__((optimize("O0")));
void while_wait() {
    while (*(yhccl_shm_addr + 64 * yhccl_new_rank) == 1)
        ;
}
void set_my_fin_to1() __attribute__((optimize("O0")));
void set_my_fin_to1() {
    *(yhccl_shm_addr + 64 * yhccl_new_rank) = 1;
}
int is_children_all_1() __attribute__((optimize("O0")));
int is_children_all_1() {
    if (yhccl_new_rank == 0) {
        for (int i = 1; i < yhccl_new_size / 2; i++) {
            if (*(yhccl_shm_addr + 64 * i) == 0)
                return 0;
        }
        return 1;
    } else {
        for (int i = yhccl_new_size / 2 + 1; i < yhccl_new_size; i++) {
            if (*(yhccl_shm_addr + 64 * i) == 0)
                return 0;
        }
        return 1;
    }
}

void set_all_0() __attribute__((optimize("O0")));
void set_all_0() {
    for (int i = 0; i < yhccl_new_size; i++) {
        *(yhccl_shm_addr + 64 * i) = 0;
    }
}

int is_all_1() __attribute__((optimize("O0")));
int is_all_1() {
    for (int i = 0; i < yhccl_new_size; i++) {
        if (*(yhccl_shm_addr + 64 * i) == 0) {
            return 0;
        }
    }
    return 1;
}

int rank0_find_childleader_is_1() __attribute__((optimize("O0")));
int rank0_find_childleader_is_1() {
    return *(yhccl_shm_addr + 64 * (yhccl_new_size / 2));
}

int is_init_intra = 0;
MPI_Comm yhccl_inter_comm;  // 注意只同步 rank0 所在的通信子

void barrier(int rank, int size) {
    if (is_init_intra == 0) {
        get_leader_and_shm_init(rank, size);

        // leader和非leader通信子划分
        char* color = (char*)calloc(size, sizeof(char));
        int i = 0;
        for (; i < leader_count; i++) {
            color[leader_rank[i]] = 1;
        }

        PMPI_Comm_split(MPI_COMM_WORLD, color[rank], rank, &yhccl_inter_comm);

        is_init_intra = 1;
    }

    if (yhccl_new_size < 16) {
        // 单层

        // 节点内
        set_my_fin_to1();

        if (yhccl_new_rank == 0) {
            while (is_all_1() == 0)
                ;
        }

        // 节点间
        if (yhccl_new_rank == 0) {
            glex_inter_barrier(yhccl_inter_comm);
            set_all_0();
        }

        // 节点内
        while_wait();
    } else {
        // 双层 2 leader
        if (yhccl_new_rank == 0 || yhccl_new_rank == yhccl_new_size / 2) {
            while (is_children_all_1() == 0)
                ;

            if (yhccl_new_rank == 0) {
                while (rank0_find_childleader_is_1() == 0)
                    ;
            } else {
                set_my_fin_to1();
                while_wait();
            }
        } else {
            set_my_fin_to1();
            while_wait();
        }

        // 只有rank0 leader会出来
        // 节点间
        if (yhccl_new_rank == 0) {
            glex_inter_barrier(yhccl_inter_comm);
            set_all_0();
        }
    }
}

int MPI_Init(int* argc, char*** argv) {
    return PMPI_Init(argc, argv);
}

int MPI_Barrier(MPI_Comm comm) {
    int rank, size;
    PMPI_Comm_rank(comm, &rank);
    PMPI_Comm_size(comm, &size);
    int gsize;
    PMPI_Comm_size(MPI_COMM_WORLD, &gsize);

    if (gsize == size)
        barrier(rank, size);
    else
        PMPI_Barrier(comm);
    return 0;
}

int MPI_Finalize(void) {
    glex_inter_destroy_ep();
    destory_shm_fd();
    return PMPI_Finalize();
}