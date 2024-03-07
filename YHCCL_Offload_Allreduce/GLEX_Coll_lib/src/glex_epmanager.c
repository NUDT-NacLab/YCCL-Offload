#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "glexcoll.h"

int leaderN = 4;
int leader_step = 4;
int NUMA_Aware_leader = 0;
int am_i_leader()
{
    if (NUMA_Aware_leader == 1)
    {
        if (intra_rank % leader_step == 0 && intra_rank / leader_step < leaderN)
        {
            return 1;
        }
    }
    if (NUMA_Aware_leader == 0)
    {
        if (intra_rank < leaderN)
            return 1;
    }
    return 0;
}
int get_my_leader_id()
{
    // printf("NUMA_Aware_leader=%d\n",NUMA_Aware_leader);
    if (NUMA_Aware_leader)
    {
        if (intra_rank % leader_step == 0 && intra_rank / leader_step < leaderN)
        {
            return intra_rank / leader_step;
        }
    }
    if (NUMA_Aware_leader == 0)
    {
        if (intra_rank < leaderN)
            return intra_rank;
    }
    puts("错误的进程想要获取leaderid");
    exit(0);
    return -1;
}
int get_ith_leader(int i)
{
    if (NUMA_Aware_leader)
    {
        return i * leader_step;
    }
    if (NUMA_Aware_leader == 0)
    {
        if (i >= leaderN)
        {
            puts("识图越界获取leader");
        }
        return i;
    }
}
void GLEX_ep_init()
{

    char *type = getenv("NUMA_Aware_leaders");
    NUMA_Aware_leader = -1;
    if (strcmp(type, "YES") == 0)
        NUMA_Aware_leader = 1;
    else if (strcmp(type, "NO") == 0)
        NUMA_Aware_leader = 0;
    else
    {
        puts("error 62");
    }

    // if (intra_rank == 0 && inter_procn > 1)
    //  printf("global_rank = %d inter_procn=%d intra_rank=%d\n",global_rank, inter_procn,intra_rank);
#ifdef PJT_NEW_VERSION
    if (inter_procn > 1 && am_i_leader())
#else
    // if (inter_procn > 1 && intra_rank % 4 == 0 && intra_rank / 4 < 6)
    if (inter_procn > 1 && intra_rank == 0)
#endif
    {
        // printf("inter_procn=%d intra_rank=%d\n", inter_procn, intra_rank);
        //每个节点上只有k个端口
        _GLEXCOLL.ep_attr.type = GLEX_EP_TYPE_FAST;
        _GLEXCOLL.ep_attr.mpq_type = GLEX_MPQ_TYPE_HIGH_CAPACITY;
        _GLEXCOLL.ep_attr.eq_type = GLEX_EQ_TYPE_HIGH_CAPACITY;
        _GLEXCOLL.ep_attr.key = 1;
        _GLEXCOLL.ep_attr.dq_capacity = inter_procn * 2048;
        _GLEXCOLL.ep_attr.mpq_capacity = inter_procn * 2048;
        _GLEXCOLL.ep_attr.eq_capacity = inter_procn * 2048;
        _GLEXCOLL.ep_attr.num = 1;
        // get_my_leader_id(); //_GLEXCOLL.intra_rank / CorePerNuma;

        // printf("grank = %d ep num= %d\n",_GLEXCOLL.global_rank,_GLEXCOLL.ep_attr.num);
        //创建软件端点

        int ret = glex_create_ep(_GLEXCOLL.dev, &(_GLEXCOLL.ep_attr), &(_GLEXCOLL.ep));
        if (ret != GLEX_SUCCESS)
        {
            printf("_create_ep() get_my_leader_id()=%d, return: %d\n", get_my_leader_id(), ret);
            exit(1);
        }
        _GLEXCOLL.ep_credit = EP_CREDIT_MAX;
        _GLEXCOLL.event_credit = Event_CREDIT_MAX;
        // AllGather获取每个leader的端点地址
        glex_ep_addr_t my_ep_addr;
        _GLEXCOLL.ep_addrs = (glex_ep_addr_t *)malloc(sizeof(*(_GLEXCOLL.ep_addrs)) * inter_procn);
        glex_get_ep_addr(_GLEXCOLL.ep, &my_ep_addr);
        MPI_Allgather((void *)&my_ep_addr, sizeof(my_ep_addr), MPI_CHAR, _GLEXCOLL.ep_addrs, sizeof(my_ep_addr), MPI_CHAR, Comm_inter);
        // printf("grank = %d, my_ep_addr = %#llx\n",_GLEXCOLL.global_rank,my_ep_addr.v);
        // printf("root nid_id = 0x%x,root ep_addr = %#llx\n \n",_GLEXCOLL.dev_attr.nic_id,(long long)my_ep_addr.v);
        // for(int i = 0;i<_GLEXCOLL.inter_procn;i++)
        // {
        // 	printf("%#llx \t",_GLEXCOLL.ep_addrs[i].v);
        // }
        // puts("");

        //设置allreduce root进程，其intra_rank ==0 且inter_rank == 0
        if (inter_rank == 0)
        {
            //根节点leader进程
            _GLEXCOLL._allreduce.root = inter_rank;
            MPI_Bcast(&(_GLEXCOLL._allreduce.root), 1, MPI_INT, 0, Comm_inter);
        }
        else
        {
            //叶节点leader进程
            MPI_Bcast(&(_GLEXCOLL._allreduce.root), 1, MPI_INT, 0, Comm_inter);
        }
        // printf("grank = %d,inter_rank = %d,intra_rank = %d,root = %d\n",_GLEXCOLL.global_rank,_GLEXCOLL.inter_rank,_GLEXCOLL.intra_rank,_GLEXCOLL.allreduce_root);

        // for(int i = 0;i<_GLEXCOLL.inter_procn;i++)
        // {
        // 	printf("ep_addr %d = %#llx \n",i,_GLEXCOLL.ep_addrs[i].v);
        // }
        //做GLEX规约请求的部分初始化
        GLEX_Coll_req.next = NULL;
        GLEX_Coll_req.data = GLEXCOLL_databuf;
        GLEX_Coll_req.len = 64;
        GLEX_Coll_req_backup.next = NULL;
        GLEX_Coll_req_backup.data = GLEXCOLL_databuf + 128;
        GLEX_Coll_req_backup.len = 64;
    }

    // if (inter_procn > 1 && intra_rank == 0)
    // {
    //     //每节点只有一个进程，0号进程需要开启两个端口
    //     _GLEXCOLL.ep_attr1.type = GLEX_EP_TYPE_FAST;
    //     _GLEXCOLL.ep_attr1.mpq_type = GLEX_MPQ_TYPE_HIGH_CAPACITY;
    //     _GLEXCOLL.ep_attr1.eq_type = GLEX_EQ_TYPE_HIGH_CAPACITY;
    //     _GLEXCOLL.ep_attr1.key = 1 + intra_rank;
    //     _GLEXCOLL.ep_attr1.dq_capacity = inter_procn * 1024;
    //     _GLEXCOLL.ep_attr1.mpq_capacity = inter_procn * 1024;
    //     _GLEXCOLL.ep_attr1.eq_capacity = inter_procn * 1024;
    //     _GLEXCOLL.ep_attr1.num = 6;

    //     int ret = glex_create_ep(_GLEXCOLL.dev, &(_GLEXCOLL.ep_attr1), &(_GLEXCOLL.ep1));
    //     if (ret != GLEX_SUCCESS)
    //     {
    //         printf("_create_ep(), return: %d\n", ret);
    //         exit(1);
    //     }

    //     glex_ep_addr_t my_ep_addr;
    //     _GLEXCOLL.ep_addrs1 = (glex_ep_addr_t *)malloc(sizeof(*(_GLEXCOLL.ep_addrs1)) * inter_procn);
    //     glex_get_ep_addr(_GLEXCOLL.ep1, &my_ep_addr);
    //     // MPI_Allgather((void *)&my_ep_addr, sizeof(my_ep_addr), MPI_CHAR, _GLEXCOLL.ep_addrs1, sizeof(my_ep_addr), MPI_CHAR, Comm_inter);
    // }
    // puts("check finish ep init");
}
void GLEX_ep_destroy()
{
#ifdef PJT_NEW_VERSION
    if (inter_procn > 1 && am_i_leader())
#else
    // if (inter_procn > 1 && intra_rank % 4 == 0 && intra_rank / 4 < 6)
#endif
        if (inter_procn > 1 && intra_rank == 0)
        {
            //销毁通信端口
            if (inter_procn > 1)
            {
                glex_ret_t ret;
                ret = glex_destroy_ep(_GLEXCOLL.ep);
                if (ret != GLEX_SUCCESS)
                {
                    printf("_destroy_ep(), return: %d\n", ret);
                    exit(1);
                }
                // if (inter_procn > 1 && intra_rank == 0)
                // {
                //     ret = glex_destroy_ep(_GLEXCOLL.ep1);
                //     // puts("check");
                //     if (ret != GLEX_SUCCESS)
                //     {
                //         printf("_destroy_ep(), return: %d\n", ret);
                //         exit(1);
                //     }
                // }
                ret = glex_close_device(_GLEXCOLL.dev);
                if (ret != GLEX_SUCCESS)
                {
                    printf("_close_device(), return: %d\n", ret);
                    exit(1);
                }
                free(_GLEXCOLL.ep_addrs);
            }
        }
}