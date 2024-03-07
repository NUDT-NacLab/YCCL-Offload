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
#include <math.h>
#include "glexcoll.h"
#include "BasicDataStructure.h"
#include "glexallreduce.h"

int K_nominal_tree_stepN(int procn, int k);
int K_nominal_tree_my_stepN(int rank, int procn, int k);

int mmin(int a, int b)
{
	return a < b ? a : b;
}
//////////////////////构造规约树/////////////////////////////////////////////////////////
//规约树的信息
struct TreeLevelInfo _ReduceTreeS[20];
int TreeLevelN;
int TreeNumber;
char host_name[MPI_MAX_PROCESSOR_NAME] = "PTT";

union glex_coll_data_loc data_loc[GLEX_COLL_REDUCE_MAX_DATA_UNITS];

struct glex_imm_mp_req GLEX_Coll_req, GLEX_Coll_req_backup;
char GLEXCOLL_databuf[1024];
uint32_t GLEXCOLL_databuf_len;

union _GLEXCOLL_INIT_K_BINOMIALT_data
{
	/* data */
	int data[10];
	struct _tmp
	{
		/* data */
		int father;
		int start;
		int end;
		int childn;
		int childRank;
		int topology_level;
	} T;
} tree_data_trans, tree_data_transTmp, tree_data_transV[500];
static int tree_data_transN;
void PrintTreeLevelInfo()
{
	//现在开始输出检查每个节点的信息
	int r = 0;
	int j;
	for (j = 0; j < TreeNumber; j++)
	{
		if (inter_rank == 0)
			printf("Tree id = %d \n", j);
		MPI_Barrier(Comm_inter);
		MPI_Barrier(Comm_inter);
		for (r = 0; r < inter_procn; r++)
		{
			if (inter_rank == r)
			{
				//输出自身节点的层次信息
				printf("inter_rank = %d  my_ep_addr =  %#llx\n", inter_rank, _GLEXCOLL.ep_addrs[inter_rank].v);
				int i;
				switch (_ReduceTreeS[j].type)
				{
				case LEAF:
					printf("LEAF\t:brothersN = %d,DownReduceRank = %d\n",
						   _ReduceTreeS[j].brothersN, _ReduceTreeS[j].DownReduceRank);
					printf("parent addr =  %#llx\n", _ReduceTreeS[j].parentAddr.v);
					break;
				case MID:
					printf("MID\t:brothersN = %d,DownReduceRank = %d childsN =%d\n",
						   _ReduceTreeS[j].brothersN, _ReduceTreeS[j].DownReduceRank, _ReduceTreeS[j].childsN);
					printf("parent addr =  %#llx\n", _ReduceTreeS[j].parentAddr.v);
					for (i = 0; i < _ReduceTreeS[j].childsN; i++)
					{
						printf("type:%d %#llx\t", _ReduceTreeS[j].childTypes[i], _ReduceTreeS[j].childAddrs[i].v);
					}
					puts("");
					break;
				case ROOT:
					printf("ROOT:\t DownReduceRank = %d,childsN =%d\n", _ReduceTreeS[j].DownReduceRank, _ReduceTreeS[j].childsN);
					for (i = 0; i < _ReduceTreeS[j].childsN; i++)
					{
						printf("type:%d %#llx\t", _ReduceTreeS[j].childTypes[i], _ReduceTreeS[j].childAddrs[i].v);
					}
					puts("");
					break;
				default:
					break;

					puts("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
				}
				puts("-----------------------------------------------------------------------------------------------------------------------------------------");

				fflush(stdout);
			}
			MPI_Barrier(Comm_inter);
			MPI_Barrier(Comm_inter);
		}
	}
}
int stringCmp(const void *a, const void *b)
{
	return strcmp(a, b);
}

//对节点进行重标号的初始化方法
int P2PDistance(int source, int dest)
{
	//最简单的距离就是逻辑距离
	return abs(source - dest);
}
int mmax(int a, int b)
{
	return a > b ? a : b;
}
int TwoAryReduceDist(int parent, int child1, int child2)
{
	return mmax(P2PDistance(parent, child1), P2PDistance(parent, child2));
}
int search_minimum_between(int start, int end)
{
	if (end < start)
		return -1;
	if (end == start)
		return start;
	int max = 0;
	int mid = 0;
	int startB = start;
	int endB = end;
	while (end > start)
	{
		mid = ((end + start) / 2);
		// printf("%d\n",mid);
		if (mid + 1 <= end)
		{
			if (TwoAryReduceDist(mid, startB, endB) < TwoAryReduceDist(mid + 1, startB, endB))
			{
				//极小值点在[start,mid]
				// puts("check");
				end = mid;
			}
			else
			{
				start = mid + 1;
			}
		}
		else
		{
		}
	}
	return end;
}

int get_child_N(int start, int end, int step)
{
	int K = 2;
	end = mmin(end, start + step);
	if (end - start > K)
		return K;
	else
		return end - start;
}

int Get_Ith_child(int start, int end, int step, int i)
{
	int n = end - start;
	int s = n / 2;
	int remain = n % s;
	if (i == 0)
		return start + 1;
	else if (i == 1)
		return mmin(start + step, start + 1 + s + remain);
	else
	{
		puts("error child");
	}
}

int Get_Ith_child_type(int start, int end, int step, int i)
{
	int n = end - start;
	int s = n / 2;
	int remain = n % s;
	if (i == 0)
	{
		//此时是左子树。左子树的范围是[start+1,Get_Ith_child(1)]
		if (start + 1 == Get_Ith_child(start, end, step, 1) - 1)
			return LEAF;
		else
			return MID;

		return start + 1;
	}
	else if (i == 1)
	{
		//此时是右子树，右子树的范围是[Get_Ith_child(1),end]
		if (Get_Ith_child(start, end, step, 1) == end)
			return LEAF;
		else
			return MID;
	}
	else
	{
		puts("error child");
		exit(0);
	}
}

int ChildNumber(int start, int end, int K)
{
	if (end - start + 1 <= K)
		return (end - start + 1);
	return K;
}
int IthChild(int start, int end, int K, int i)
{
	if (end - start + 1 <= K)
		return start + i;
	int step = (end - start + 1) / K;
	int remain = (end - start + 1) % K;
	return start + i * step + (i < remain ? i : remain);
}
void _GLEXCOLL_Init_ReduceTree_Rack_aware()
{
	//此处的核心目的是构建机架感知的规约树
	//第一步读取文件
	// puts("check");
	char filename[200] = "PJT_Partition_adviseFile.txt";
	// system("pwd");
	FILE *fp = fopen(filename, "r");
	if (fp == 0)
		puts("error file  PJT_Partition_adviseFile");
	int len;
	fscanf(fp, "%d", &len);
	int *partition_vec = (int *)malloc(sizeof(int) * len);
	for (int i = 0; i < len; ++i)
	{
		fscanf(fp, "%d", &(partition_vec[i]));
		// if(global_rank == 0)
		// 	printf("%d %d\n",len,partition_vec[i]);
	}
	fclose(fp);
	// if(global_rank == 0)
	// 	for(int i = 0;i<inter_procn;i++)
	// 		printf("%d\n",host_ids[i]);
	//根据host_ids初始化虚拟机架的开始inter_rank号
	int Rack_PositionStart[1000]; // ranck[进程顺序机架号] = rack开始进程inter_rank
	int Rack_PositionEnd[1000];	  // ranck[进程顺序机架号] = rack结束进程inter_rank
	int *Rank_To_Rackid = (int *)malloc(sizeof(int) * inter_procn);
	int Rack_Num = 0;
	int rackend = 0;
	for (int i = 0; i < inter_procn; i++)
	{
		if (host_ids[i] >= partition_vec[rackend])
		{
			//跨越了一个新的机架
			Rack_PositionStart[Rack_Num++] = i;
			while (host_ids[i] >= partition_vec[rackend])
				rackend++;
		}
	}
	for (int i = 0; i < Rack_Num - 1; i++)
	{
		Rack_PositionEnd[i] = Rack_PositionStart[i + 1] - 1;
	}
	Rack_PositionEnd[Rack_Num - 1] = inter_procn - 1;

	for (int i = 0; i < Rack_Num; i++)
	{
		// if(global_rank == 0)
		// 	printf("(%d,%d)\n",Rack_PositionStart[i],Rack_PositionEnd[i]);
		for (int a = Rack_PositionStart[i]; a <= Rack_PositionEnd[i]; a++)
		{
			Rank_To_Rackid[a] = i;
		}
	}
	// if(global_rank == 0)
	// 	for(int i = 0;i< Rack_Num;i++)
	// 	{
	// 		printf("inter_rank=%d host_ids=%d\n", Rack_PositionStart[i],host_ids[Rack_PositionStart[i]]);
	// 	}
	int K_inter = 2;
	int K_intra = 2;
	//开始划分
	int root = 0;
	TreeNumber = 0;
	if (inter_rank == root)
	{
		int start = 1;
		int end = inter_procn - 1;
		_ReduceTreeS[TreeNumber].type = ROOT;
		_ReduceTreeS[TreeNumber].parentAddr = _GLEXCOLL.ep_addrs[inter_rank];
		_ReduceTreeS[TreeNumber].parentID = root;
		_ReduceTreeS[TreeNumber].brothersN = 0;
		_ReduceTreeS[TreeNumber].DownReduceRank = 0;
		tree_data_transTmp.T.father = inter_rank;
		if (Rack_Num > 1)
		{
			//机架间
			int childN = ChildNumber(1, Rack_Num - 1, K_inter);
			//先看下本地机架节点数量是否大于1
			int shiftChildN = 0;
			// printf("%d %d\n",Rack_PositionEnd[0],Rack_PositionStart[0]);
			if (Rack_PositionEnd[0] - Rack_PositionStart[0] >= 1)
			{
				//本地机架至少两个节点
				shiftChildN = 1;
				tree_data_transTmp.T.childn = childN + shiftChildN;
				tree_data_transTmp.T.childRank = 1;
				tree_data_transTmp.T.start = Rack_PositionStart[0] + 1;
				tree_data_transTmp.T.end = Rack_PositionEnd[0];

				_ReduceTreeS[TreeNumber].childAddrs[0] = _GLEXCOLL.ep_addrs[tree_data_transTmp.T.start];
				_ReduceTreeS[TreeNumber].childIds[0] = tree_data_transTmp.T.start;
				// _ReduceTreeS[TreeNumber].childAddrs1[0] = _GLEXCOLL.ep_addrs1[tree_data_transTmp.T.start];
				if (tree_data_transTmp.T.start == tree_data_transTmp.T.end)
					_ReduceTreeS[TreeNumber].childTypes[0] = LEAF;
				else
					_ReduceTreeS[TreeNumber].childTypes[0] = MID;
				MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, tree_data_transTmp.T.start, 996, Comm_inter);
				// printf("childn=%d childRank=%d start=%d end=%d inter_rank=%d\n",
				// 			tree_data_transTmp.T.childn,
				// 			tree_data_transTmp.T.childRank,
				// 			tree_data_transTmp.T.start,
				// 			tree_data_transTmp.T.end,
				// 			inter_rank);
			}
			_ReduceTreeS[TreeNumber].childsN = childN + shiftChildN;
			for (int i = 0; i < childN; i++)
			{
				int child_start = Rack_PositionStart[IthChild(1, Rack_Num - 1, K_inter, i)];
				int child_end = inter_procn - 1;
				if (i != childN - 1)
					child_end = Rack_PositionStart[IthChild(1, Rack_Num - 1, K_inter, i + 1)] - 1;
				_ReduceTreeS[TreeNumber].childAddrs[i + shiftChildN] = _GLEXCOLL.ep_addrs[child_start];
				// _ReduceTreeS[TreeNumber].childAddrs1[i + shiftChildN] = _GLEXCOLL.ep_addrs1[child_start];
				_ReduceTreeS[TreeNumber].childIds[i + shiftChildN] = child_start;
				if (child_start == child_end)
					_ReduceTreeS[TreeNumber].childTypes[i + shiftChildN] = LEAF;
				else
					_ReduceTreeS[TreeNumber].childTypes[i + shiftChildN] = MID;
				tree_data_transTmp.T.childn = childN + shiftChildN;
				tree_data_transTmp.T.childRank = 1 + i + shiftChildN;
				tree_data_transTmp.T.start = child_start;
				tree_data_transTmp.T.end = child_end;
				MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, tree_data_transTmp.T.start, 996, Comm_inter);

				// printf("childn=%d childRank=%d start=%d end=%d inter_rank=%d\n",
				// 			tree_data_transTmp.T.childn,
				// 			tree_data_transTmp.T.childRank,
				// 			tree_data_transTmp.T.start,
				// 			tree_data_transTmp.T.end,
				// 			inter_rank);
			}
		}
		else
		{
			//机架内
			int start = 1;
			int end = inter_procn - 1;
			int childN = ChildNumber(start, end, K_intra);
			_ReduceTreeS[TreeNumber].childsN = childN;

			for (int i = 0; i < childN; i++)
			{
				int child_start = IthChild(start, end, K_intra, i);
				int child_end = end;
				if (i != childN - 1)
					child_end = IthChild(start, end, K_intra, i + 1) - 1;
				_ReduceTreeS[TreeNumber].childAddrs[i] = _GLEXCOLL.ep_addrs[child_start];
				// _ReduceTreeS[TreeNumber].childAddrs1[i] = _GLEXCOLL.ep_addrs1[child_start];
				_ReduceTreeS[TreeNumber].childIds[i] = child_start;
				if (child_start == child_end)
					_ReduceTreeS[TreeNumber].childTypes[i] = LEAF;
				else
					_ReduceTreeS[TreeNumber].childTypes[i] = MID;
				tree_data_transTmp.T.start = child_start;
				tree_data_transTmp.T.end = child_end;
				tree_data_transTmp.T.childn = childN;
				tree_data_transTmp.T.childRank = 1 + i;
				MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, tree_data_transTmp.T.start, 996, Comm_inter);

				// printf("childn=%d childRank=%d start=%d end=%d inter_rank=%d\n",
				// 			tree_data_transTmp.T.childn,
				// 			tree_data_transTmp.T.childRank,
				// 			tree_data_transTmp.T.start,
				// 			tree_data_transTmp.T.end,
				// 			inter_rank);
			}
		}
	}
	else
	{
		MPI_Status status;
		//接收来自父节点的数据分组
		// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		MPI_Recv(&tree_data_trans, sizeof(tree_data_trans), MPI_CHAR, MPI_ANY_SOURCE, 996, Comm_inter, &status);
		// printf(" parent_rank = %d brotherN=%d UpReduceRank=%d start=%d end=%d inter_rank=%d rackid = %d\n",
		// 			tree_data_trans.T.father,
		// 			tree_data_trans.T.childn,
		// 			tree_data_trans.T.childRank,
		// 			tree_data_trans.T.start,
		// 			tree_data_trans.T.end,
		// 			inter_rank,
		// 			Rank_To_Rackid[inter_rank]);
		int startrank = tree_data_trans.T.start;
		int endrank = tree_data_trans.T.end;
		if (startrank == endrank)
		{
			//叶子节点，中止
			_ReduceTreeS[TreeNumber].type = LEAF;
		}
		else
		{
			_ReduceTreeS[TreeNumber].type = MID;
		}
		_ReduceTreeS[TreeNumber].parentAddr = _GLEXCOLL.ep_addrs[tree_data_trans.T.father];
		_ReduceTreeS[TreeNumber].parentID = tree_data_trans.T.father;
		_ReduceTreeS[TreeNumber].brothersN = tree_data_trans.T.childn;
		_ReduceTreeS[TreeNumber].DownReduceRank = tree_data_trans.T.childRank;
		if (startrank != endrank)
		{
			tree_data_transTmp.T.father = inter_rank;
			//只处理非叶节点
			//接下来处理孩子的数量和孩子节点的类型
			int rack_start = Rank_To_Rackid[startrank];
			int rack_end = Rank_To_Rackid[endrank];
			// printf("rack_start = %d,rack_end = %d\n",rack_start,rack_end);
			if (rack_start != rack_end)
			{
				//机架间
				//第一步计算孩子的数量
				int childN = ChildNumber(rack_start + 1, rack_end, K_inter);
				int shiftChildN = 0;
				if (Rack_PositionEnd[rack_start] - startrank >= 1)
				{
					//本地机架至少两个节点
					shiftChildN = 1;
					tree_data_transTmp.T.childn = childN + shiftChildN;
					tree_data_transTmp.T.childRank = 1;
					tree_data_transTmp.T.start = startrank + 1;
					tree_data_transTmp.T.end = Rack_PositionEnd[rack_start];
					_ReduceTreeS[TreeNumber].childAddrs[0] = _GLEXCOLL.ep_addrs[tree_data_transTmp.T.start];
					// _ReduceTreeS[TreeNumber].childAddrs1[0] = _GLEXCOLL.ep_addrs1[tree_data_transTmp.T.start];
					_ReduceTreeS[TreeNumber].childIds[0] = tree_data_transTmp.T.start;
					if (tree_data_transTmp.T.start == tree_data_transTmp.T.end)
						_ReduceTreeS[TreeNumber].childTypes[0] = LEAF;
					else
						_ReduceTreeS[TreeNumber].childTypes[0] = MID;

					MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, tree_data_transTmp.T.start, 996, Comm_inter);
					// printf("BBB childn=%d childRank=%d start=%d end=%d inter_rank=%d\n",
					// 			tree_data_transTmp.T.childn,
					// 			tree_data_transTmp.T.childRank,
					// 			tree_data_transTmp.T.start,
					// 			tree_data_transTmp.T.end,
					// 			inter_rank);
				}
				_ReduceTreeS[TreeNumber].childsN = childN + shiftChildN;
				//接下来处理机架间，将机架间的节点进行分割
				for (int i = 0; i < childN; i++)
				{
					int child_start = Rack_PositionStart[IthChild(rack_start + 1, rack_end, K_inter, i)];
					int child_end = endrank;
					if (i != childN - 1)
						child_end = Rack_PositionStart[IthChild(rack_start + 1, rack_end, K_inter, i + 1)] - 1;
					_ReduceTreeS[TreeNumber].childAddrs[i + shiftChildN] = _GLEXCOLL.ep_addrs[child_start];
					// _ReduceTreeS[TreeNumber].childAddrs1[i + shiftChildN] = _GLEXCOLL.ep_addrs1[child_start];
					_ReduceTreeS[TreeNumber].childIds[i + shiftChildN] = child_start;
					if (child_start == child_end)
						_ReduceTreeS[TreeNumber].childTypes[i + shiftChildN] = LEAF;
					else
						_ReduceTreeS[TreeNumber].childTypes[i + shiftChildN] = MID;
					tree_data_transTmp.T.childn = childN + shiftChildN;
					tree_data_transTmp.T.childRank = 1 + i + shiftChildN;
					tree_data_transTmp.T.start = child_start;
					tree_data_transTmp.T.end = child_end;
					MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, tree_data_transTmp.T.start, 996, Comm_inter);

					// printf("BBB childn=%d childRank=%d start=%d end=%d inter_rank=%d\n",
					// 			tree_data_transTmp.T.childn,
					// 			tree_data_transTmp.T.childRank,
					// 			tree_data_transTmp.T.start,
					// 			tree_data_transTmp.T.end,
					// 			inter_rank);
				}
			}
			else
			{
				//机架内
				// printf("chekc 机架内 rank = %d\n",inter_rank);
				int start = inter_rank + 1;
				int end = endrank;
				int childN = ChildNumber(start, end, K_intra);
				_ReduceTreeS[TreeNumber].childsN = childN;

				for (int i = 0; i < childN; i++)
				{
					int child_start = IthChild(start, end, K_intra, i);
					int child_end = endrank;
					if (i != childN - 1)
						child_end = IthChild(start, end, K_intra, i + 1) - 1;
					_ReduceTreeS[TreeNumber].childAddrs[i] = _GLEXCOLL.ep_addrs[child_start];
					// _ReduceTreeS[TreeNumber].childAddrs1[i] = _GLEXCOLL.ep_addrs1[child_start];
					_ReduceTreeS[TreeNumber].childIds[i] = child_start;
					if (child_start == child_end)
						_ReduceTreeS[TreeNumber].childTypes[i] = LEAF;
					else
						_ReduceTreeS[TreeNumber].childTypes[i] = MID;
					tree_data_transTmp.T.start = child_start;
					tree_data_transTmp.T.end = child_end;
					tree_data_transTmp.T.childn = childN;
					tree_data_transTmp.T.childRank = 1 + i;
					MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, tree_data_transTmp.T.start, 996, Comm_inter);

					// printf("childn=%d childRank=%d start=%d end=%d inter_rank=%d\n",
					// 			tree_data_transTmp.T.childn,
					// 			tree_data_transTmp.T.childRank,
					// 			tree_data_transTmp.T.start,
					// 			tree_data_transTmp.T.end,
					// 			inter_rank);
				}
			}
		}
	}
	// for(int r = 1;r<inter_procn;r++)
	// {
	// 	if(inter_rank == r)
	// 	{
	// 		printf(" parent_rank = %d brotherN=%d UpReduceRank=%d start=%d end=%d inter_rank=%d rackid = %d\n",
	// 			   tree_data_trans.T.father,
	// 			   tree_data_trans.T.childn,
	// 			   tree_data_trans.T.childRank,
	// 			   tree_data_trans.T.start,
	// 			   tree_data_trans.T.end,
	// 			   inter_rank,
	// 			   Rank_To_Rackid[inter_rank]);
	// 	}
	// 		MPI_Barrier(Comm_inter);
	// }
	free(partition_vec);
	TreeLevelN = 1;
	TreeNumber = 1;
	// PrintTreeLevelInfo();

	// exit(0);
}

void _GLEXCOLL_Init_TJ3_topology_aware_ary_tree()
{
	int get_childN(int start, int end, int *hierarchyVec, int hierarchyn, int *childv)
	{
		// start是parent [start,end)
		//  end = mmin(inter_procn, end);
		int childn = 0;
		int child = start + 1;
		int breakflag = 0;
		int step = 1;
		for (int h = 0; h < hierarchyn; h++)
		{
			step *= hierarchyVec[h];
			//对每个层次的孩子节点
			for (int t = 1; t < hierarchyVec[h + 1]; t++)
			{

				if (child < end)
				{
					childv[childn] = child;
					// printf("%d child %d = %d\n", inter_rank, childn, childv[childn]);
					childn++;
				}
				else
					return childn;
				child += step;
			}
		}
	}

	int root = 0;
	TreeNumber = 0;
	TreeLevelN = 0;
	//先初始化一下_ReduceTreeS
	for (int i = 0; i < 2; i++)		  //最大20层树高
		for (int j = 0; j < 256; j++) //每层最多256个child
			_ReduceTreeS[i].childTypes[j] = MID;
	_ReduceTreeS[0].childsN = 0;   //当前节点的孩子节点数量
	_ReduceTreeS[0].brothersN = 0; //当前节点的兄弟节点数量
	int hierarchyN = 30;
	int hierarchyVec[hierarchyN + 1];
	hierarchyVec[0] = 1;
	hierarchyVec[1] = 8;
	hierarchyVec[2] = 16;
	hierarchyVec[3] = 2;
	for (int i = 4; i < hierarchyN + 1; i++)
		hierarchyVec[i] = 6;

	if (inter_rank == 0)
	{
		// root
		_ReduceTreeS[TreeLevelN].type = ROOT;
		//记录下父节点端口地址信息,也就是其自身。
		_ReduceTreeS[TreeLevelN].parentAddr = _GLEXCOLL.ep_addrs[root];
		_ReduceTreeS[TreeLevelN].parentID = root;
		int dist = 1;
		int start = 0;
		int childn = get_childN(0, inter_procn, hierarchyVec, hierarchyN, _ReduceTreeS[TreeLevelN].childIds);
		_ReduceTreeS[TreeLevelN].childsN = childn;
		_ReduceTreeS[TreeLevelN].DownReduceRank = 0;
		for (int i = 0; i < childn; i++)
		{
			int child = _ReduceTreeS[TreeLevelN].childIds[i];
			_ReduceTreeS[TreeLevelN].childAddrs[i] = _GLEXCOLL.ep_addrs[child];
			// _ReduceTreeS[TreeLevelN].childAddrs1[i] = _GLEXCOLL.ep_addrs1[child];
			tree_data_transTmp.T.father = 0;
			tree_data_transTmp.T.start = child;
			if (i == childn - 1)
				tree_data_transTmp.T.end = inter_procn;
			else
				tree_data_transTmp.T.end = _ReduceTreeS[TreeLevelN].childIds[i + 1];
			tree_data_transTmp.T.childn = childn;
			tree_data_transTmp.T.childRank = i + 1;
			MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, child, 0, Comm_inter);
		}
	}
	else
	{
		MPI_Status s;
		MPI_Recv(&tree_data_trans, sizeof(tree_data_trans), MPI_CHAR, MPI_ANY_SOURCE, 0, Comm_inter, &s);
		int start = tree_data_trans.T.start;
		int end = tree_data_trans.T.end;
		int parent = tree_data_trans.T.father;
		// 1.初始化节点类型
		if (end == start + 1)
		{
			_ReduceTreeS[TreeLevelN].type = LEAF;
			_ReduceTreeS[TreeLevelN].childsN = 0;
		}
		else
			_ReduceTreeS[TreeLevelN].type = MID;
		// 2. 设置parent地址和rank
		_ReduceTreeS[TreeLevelN].parentID = parent;
		_ReduceTreeS[TreeLevelN].parentAddr = _GLEXCOLL.ep_addrs[parent];
		// 3. 设置兄弟数量
		_ReduceTreeS[TreeLevelN].brothersN = tree_data_trans.T.childn;
		// 4. 设置自身进程在一次规约规约时的序号
		_ReduceTreeS[TreeLevelN].DownReduceRank = tree_data_trans.T.childRank;

		if (_ReduceTreeS[TreeLevelN].type == MID)
		{
			int childn = get_childN(start, end, hierarchyVec, hierarchyN, _ReduceTreeS[TreeLevelN].childIds);
			_ReduceTreeS[TreeLevelN].childsN = childn;
			for (int i = 0; i < childn; i++)
			{
				int child = _ReduceTreeS[TreeLevelN].childIds[i];
				_ReduceTreeS[TreeLevelN].childAddrs[i] = _GLEXCOLL.ep_addrs[child];
				// _ReduceTreeS[TreeLevelN].childAddrs1[i] = _GLEXCOLL.ep_addrs1[child];
				tree_data_transTmp.T.father = inter_rank;
				tree_data_transTmp.T.start = child;
				if (i == childn - 1)
					tree_data_transTmp.T.end = end;
				else
					tree_data_transTmp.T.end = _ReduceTreeS[TreeLevelN].childIds[i + 1];
				tree_data_transTmp.T.childn = childn;
				tree_data_transTmp.T.childRank = i + 1;
				MPI_Send(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, child, 0, Comm_inter);
			}
		}
	}
	// printf("rank=%d type = %d, ")
	MPI_Barrier(Comm_inter);
	// fflush(stdout);
	// exit(0);
}
void _GLEXCOLL_Init_ReduceTree()
{

	//先初始化一下_ReduceTreeS
	for (int i = 0; i < 20; i++)
		for (int j = 0; j < 32; j++)
			_ReduceTreeS[i].childTypes[j] = MID;

	int algorithm = 2;
	char *type = getenv("INTER_ALLREDUCE_TYPE");
	if (type == 0)
	{
		puts("_GLEXCOLL_Init_ReduceTree type=0");
		exit(1);
	}
	if (strcmp(type, "ORIGINAL") == 0)
	{
		//基础节点间算法
		algorithm = 1;
	}
	else if (strcmp(type, "PERFORMANCE_AWARENESS") == 0)
	{
		//
		algorithm = 2;
	}
	else if (strcmp(type, "TJ_3_TOPOLOGY_AWARE") == 0)
	{
		algorithm = 3;
	}
	MPI_Status localstatus;
	MPI_Request localreq, localreqV[100];
	int root = 0;
	int K = Childn_K;
	int tmp = -1;

	switch (algorithm)
	{
	case 1: /* 树形规约算法K-Trees*/
		/* code */
		{
			//用MPI的角度来从上到下初始化
			int start = 0;
			int end = 0;
			tree_data_transN = 0;
			if (inter_rank == root)
			{
				int p = 0;
				int q = inter_procn - 1;
				tree_data_transTmp.data[0] = -1;
				tree_data_transTmp.data[1] = p;
				tree_data_transTmp.data[2] = q;
				tree_data_transTmp.T.childn = (q - p + 1) / K + (((q - p + 1) % K > 0) ? 1 : 0);
				tree_data_transTmp.T.childRank = -1; //这是发给规约根的消息，规约根不需要childRank,故而随意设置
				MPI_Isend(&tree_data_transTmp, sizeof(tree_data_transTmp), MPI_CHAR, inter_rank, 0, Comm_inter, &(localreqV[tree_data_transN++]));
			}
			int father;
			{
				if (MPI_Recv(&tree_data_trans, sizeof(tree_data_trans), MPI_CHAR, MPI_ANY_SOURCE, 0, Comm_inter, &localstatus) != MPI_SUCCESS)
				{
					puts("MPI RECV ERROR");
					MPI_Finalize();
					exit(0);
				}
				// printf("inter_rank = %d,recvmsg from %d\n",inter_rank,localstatus.MPI_SOURCE);

				if (tree_data_trans.T.father == -1)
				{
					//这是根进程才会受到的其中一种消息
					_ReduceTreeS[TreeLevelN].type = ROOT;
					//记录下父节点端口地址信息,也就是其自身。
					_ReduceTreeS[TreeLevelN].parentAddr = _GLEXCOLL.ep_addrs[root];
					_ReduceTreeS[TreeLevelN].parentID = root;
				}
				else
				{

					//记录下父节点端口地址信息
					_ReduceTreeS[TreeLevelN].parentAddr = _GLEXCOLL.ep_addrs[tree_data_trans.T.father];
					_ReduceTreeS[TreeLevelN].parentID = tree_data_trans.T.father;
					//记录下自己和兄弟的数量
					_ReduceTreeS[TreeLevelN].brothersN = tree_data_trans.T.childn;
					//记录下自己在上方规约中的rank位置
					_ReduceTreeS[TreeLevelN].DownReduceRank = tree_data_trans.T.childRank;
					if (tree_data_trans.T.start == tree_data_trans.T.end)
					{
						//抵达leaf
						_ReduceTreeS[TreeLevelN].type = LEAF;
						TreeLevelN++;
					}
					else
						_ReduceTreeS[TreeLevelN].type = MID;
				}
				if (_ReduceTreeS[TreeLevelN].type != LEAF)
				{
					//接下来将[start,end]区间分割为K块
					int block = 0, p, q;
					int step = (tree_data_trans.T.end - tree_data_trans.T.start) / K;
					int remain = (tree_data_trans.T.end - tree_data_trans.T.start) % K;
					// printf("myid = %d,father = %d,[%d,%d] ,step = %d,remain = %d \n",
					// inter_rank,tree_data_trans.data[0],tree_data_trans.data[1],tree_data_trans.data[2],step,remain);
					_ReduceTreeS[TreeLevelN].childsN = 0;
					int childN = (step > 0 ? K : remain);
					for (block = 0, p = tree_data_trans.T.start + 1;
						 block < K && p <= tree_data_trans.T.end;
						 block++)
					{
						//[p,q]为第block个子区域的下一步规约区间。
						q = p + (step - 1 + ((block < remain) ? 1 : 0));
						// printf("myid = %d,father = %d,[%d,%d] ->[%d,%d]\n",
						// inter_rank,tree_data_trans.data[0],tree_data_trans.T.start,tree_data_trans.T.end,p,q);
						// fflush(stdout);
						tree_data_transV[tree_data_transN].T.father = inter_rank;
						tree_data_transV[tree_data_transN].T.start = p;
						tree_data_transV[tree_data_transN].T.end = q;
						tree_data_transV[tree_data_transN].T.childn = childN;
						tree_data_transV[tree_data_transN].T.childRank = block + 1;
						int child = 0;
						//在[p,q]中选出下一步规约区间的根，
						child = p;

						// printf("inter_rank = %d [%d,%d] send to %d\n",inter_rank,p,q,child);
						// MPI_Isend(&(tree_data_transV[tree_data_transN]),sizeof(tree_data_trans),MPI_CHAR,child,0,Comm_inter,&(localreqV[tree_data_transN]));
						// tree_data_transN++;
						MPI_Send(&(tree_data_transV[tree_data_transN]), sizeof(tree_data_trans), MPI_CHAR, child, 0, Comm_inter);
						//写入叶子点相关信息
						{
							//统计叶子的数量，并记录下叶子的端口地址以方便规约
							// printf("inter rank = %d child = %d\n",inter_rank,child);
							_ReduceTreeS[TreeLevelN].childAddrs[_ReduceTreeS[TreeLevelN].childsN] = _GLEXCOLL.ep_addrs[child];
							// _ReduceTreeS[TreeLevelN].childAddrs1[_ReduceTreeS[TreeLevelN].childsN] = _GLEXCOLL.ep_addrs1[child];
							_ReduceTreeS[TreeLevelN].childIds[_ReduceTreeS[TreeLevelN].childsN] = child;
							_ReduceTreeS[TreeLevelN].childsN++;
						}
						//将这里注释掉是否会对单数据allreduce的正确性产生影响还有待评估
						// {
						// 	//找到inter_rank在下方规约中的rank号
						// 	//此处只有ROOT和MID节点一定会进入
						// 	_ReduceTreeS[TreeLevelN].DownReduceRank = 0;
						// }
						p = q + 1;
					}
					TreeLevelN++;
				}
				for (int i = 0; i < tree_data_transN; i++)
				{
					if (MPI_Wait(&(localreqV[i]), &localstatus) != MPI_SUCCESS)
					{
						puts("MPI_Wait error ");
						exit(0);
					}
				}
				// printf("inter_rank = %d,barrier finish\n",inter_rank);
			}

			MPI_Barrier(Comm_inter);
			MPI_Barrier(Comm_inter);
			if (0)
			{
				//现在开始输出检查每个节点的信息
				int r = 0;
				for (r = 0; r < inter_procn; r++)
				{
					if (inter_rank == r)
					{
						//输出自身节点的层次信息
						printf("inter_rank = %d TreeLevelN = %d my_ep_addr =  %#llx\n", inter_rank, TreeLevelN, _GLEXCOLL.ep_addrs[inter_rank].v);
						int j;
						for (j = 0; j < TreeLevelN; j++)
						{
							int i;
							switch (_ReduceTreeS[j].type)
							{
							case LEAF:
								printf("LEAF\t:brothersN = %d,DownReduceRank = %d\n",
									   _ReduceTreeS[j].brothersN, _ReduceTreeS[j].DownReduceRank);
								printf("parent addr =  %#llx\n", _ReduceTreeS[j].parentAddr.v);
								break;
							case MID:
								printf("MID\t:brothersN = %d,UpReduceRank = %d,DownReduceRank = %d childsN =%d\n",
									   _ReduceTreeS[j].brothersN, _ReduceTreeS[j].UpReduceRank, _ReduceTreeS[j].DownReduceRank, _ReduceTreeS[j].childsN);
								printf("parent addr =  %#llx\n", _ReduceTreeS[j].parentAddr.v);
								for (i = 0; i < _ReduceTreeS[j].childsN; i++)
								{
									printf("%#llx\t", _ReduceTreeS[j].childAddrs[i].v);
								}
								puts("");
								break;
							case ROOT:
								printf("ROOT:\t DownReduceRank = %d,childsN =%d\n", _ReduceTreeS[j].DownReduceRank, _ReduceTreeS[j].childsN);
								for (i = 0; i < _ReduceTreeS[j].childsN; i++)
								{
									printf("%#llx\t", _ReduceTreeS[j].childAddrs[i].v);
								}
								puts("");
								break;
							default:
								break;
							}

							puts("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
						}
						puts("-----------------------------------------------------------------------------------------------------------------------------------------");

						fflush(stdout);
					}
					MPI_Barrier(Comm_inter);
					MPI_Barrier(Comm_inter);
				}
			}
			TreeNumber = 1;
			// PrintTreeLevelInfo();
		}
		break;
	case 2: /*按照机柜延迟进行节点通信分层规约*/
		_GLEXCOLL_Init_ReduceTree_Rack_aware();
		break;
	case 3:
		/*天津2021 3号机器前两层拓扑感知规约,使用时进程数量必须为256的倍数。最佳性能为整框使用*/
		_GLEXCOLL_Init_TJ3_topology_aware_ary_tree();
	}
}
//////////////////////构造规约树结束/////////////////////////////////////////////////////////

int INTRA_ALLREDUCE_TYPE;
int ppn = 32;

extern volatile char *SHM_bufs_Recv[64];
extern volatile char *SHM_bufs_Send[64];
extern volatile char *SHM_bufs_bcast;
extern volatile char *SHM_bufs_bcast1;
void Init_SHMBUF1()
{
	int size = intra_procn * 64;
	char name[100];
	sprintf(name, "%s-%d\0", host_name, intra_rank);
	int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
	ftruncate(fd, size);
	volatile char *p = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (p == MAP_FAILED)
	{
		printf("error map\n");
		return;
	}
	for (int i = 0; i < size; ++i)
	{
		p[i] = 'S';
	}
	for (int i = 0; i < intra_procn; ++i)
	{
		SHM_bufs_Recv[i] = p + i * 64;
	}
	MPI_Barrier(Comm_intra);

	//接下来记录下每对进程之间共享内存的区域
	for (int i = 0; i < intra_procn; ++i)
	{
		sprintf(name, "%s-%d\0", host_name, i);
		fd = shm_open(name, O_RDWR, S_IRUSR | S_IWUSR);
		ftruncate(fd, size);
		SHM_bufs_Send[i] = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		SHM_bufs_Send[i] += intra_rank * 64;
		// if(intra_rank == 1)
		//     printf("%d -> %d \n",((int*)(SHM_bufs_Send[i]))[0],((int*)(SHM_bufs_Send[i]))[1] );
	}
	// printf("%c\n",*SHM_bufs_Recv[0]);
	//  *SHM_bufs_Send[0] = 'X';
	//  if(intra_rank != 0)
	//  	*SHM_bufs_Recv[0] = 'X';
	MPI_Barrier(Comm_intra);
}

void INIT_Shm_buf()
{
	int size = 64 * (intra_procn + 1);
	volatile char *p;

	//每个进程都创建一个缓冲区

	sprintf(_GLEXCOLL.name, "%s-%d\0", host_name, intra_rank);
	int fd = shm_open(_GLEXCOLL.name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
	ftruncate(fd, size);
	p = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (p == MAP_FAILED)
	{
		printf("error map\n");
		return;
	}
	memset((void *)p, 'S', size);
	MPI_Barrier(Comm_intra);
	//接下来记录下每对进程之间共享内存的区域
	for (int i = 0; i < intra_procn; ++i)
	{
		sprintf(_GLEXCOLL.name, "%s-%d\0", host_name, i);
		fd = shm_open(_GLEXCOLL.name, O_RDWR, S_IRUSR | S_IWUSR);
		ftruncate(fd, size);
		_GLEXCOLL.SHM_bufs[i] = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	}
	// if(intra_rank == 0)
	// 	for(int s = 0;s<intra_procn;s++)
	// 	{
	// 		for(int i = 0;i<intra_procn;i++)
	// 		{
	// 			char *p = _GLEXCOLL.SHM_bufs[s];
	// 			p =p +(s<<6);
	// 			printf("%c ",*p);
	// 		}
	// 		puts("");
	// 	}
}
extern void Intel_Init_SHMBUF();
void create_shemem(char *type)
{
	Init_SHMBUF1();
	// MPI_Barrier(MPI_COMM_WORLD);
	// puts("check 868");
	// printf("INTRA_ALLREDUCE_TYPE = %s\n",type);
	Intel_Init_SHMBUF();
	// MPI_Barrier(MPI_COMM_WORLD);
	// puts("check 872");
	// printf("INTRA_ALLREDUCE_TYPE = %s\n",type);
	if (strcmp(type, "REDUCE_BCAST") == 0 || strcmp(type, "REDUCE_BCAST_REORDER") == 0)
	{
		INTRA_ALLREDUCE_TYPE = 0;
		if (strcmp(type, "REDUCE_BCAST_REORDER") == 0)
			INTRA_ALLREDUCE_TYPE = 4;
		//此时所有的消息都会发送给intra leader，让leader去进行消息传输。
		int size = 64 * (intra_procn + 3);
		volatile char *p;
		// sprintf(_GLEXCOLL.name,"GLEXshm%d\0",_GLEXCOLL.global_rank/ppn);
		// printf("%s\n",_GLEXCOLL.name);
		if (intra_rank == 0)
		{
			//读取者
			int fd = shm_open(host_name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
			ftruncate(fd, size); //除了0以外63个核心。
			p = (volatile char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			if (p == MAP_FAILED)
			{
				printf("error map\n");
				return;
			}
			for (int j = 0; j < size; j++)
				p[j] = 'S';
			_GLEXCOLL.SHM_0 = _GLEXCOLL.SHM_buf = p;
			MPI_Barrier(Comm_intra);
		}
		else
		{
			//写入者
			//其它进程只能打开共享缓存文件
			MPI_Barrier(Comm_intra);
			int fd = shm_open(host_name, O_RDWR, S_IRUSR | S_IWUSR);
			ftruncate(fd, size);
			volatile char *pt = (volatile char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			if (pt == MAP_FAILED)
			{
				printf("error map\n");
				return;
			}
			p = pt + (_GLEXCOLL.intra_rank - 1) * 64;
			// printf("rank = %d %p %p\n",intra_rank,pt,p);
			_GLEXCOLL.SHM_0 = pt;
			_GLEXCOLL.SHM_buf = p;
			// printf("*p = %c\n",*(_GLEXCOLL.SHM_buf));
		}
		SHM_bufs_bcast = _GLEXCOLL.SHM_0 + (intra_procn + 1) * 64;
		SHM_bufs_bcast1 = _GLEXCOLL.SHM_0 + (intra_procn + 2) * 64;

		MPI_Barrier(Comm_intra);
	}
	else if (strcmp(type, "2_NOMINAL") == 0)
	{

		INTRA_ALLREDUCE_TYPE = 1;
		INIT_Shm_buf();

		// //puts("check");
		// INTRA_ALLREDUCE_TYPE = 1;
		// //此时所有的消息都会发送给intra leader，让leader去进行消息传输。
		// int size = 64*intra_procn;
		// volatile int *p;
		// // sprintf(_GLEXCOLL.name,"GLEXshm%d\0",_GLEXCOLL.global_rank/ppn);
		// // printf("%s\n",_GLEXCOLL.name);
		// if(intra_rank == 0)
		// {
		// 	//读取者
		// 	int fd = shm_open(host_name,O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
		// 	ftruncate(fd,size);//除了0以外63个核心。
		// 	p = (int *)mmap(NULL,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
		// 	if (p == MAP_FAILED)
		// 	{
		// 		printf("error map\n");
		// 		return ;
		// 	}
		// 	for(int i = 0;i<intra_procn;++i)
		// 	{
		// 		if(i%2 == 0)
		// 			*(p+(i<<4))= i+1;
		// 		else
		// 		{
		// 			*(p+(i<<4))= -1;
		// 		}
		// 		//printf("%d %d\n",i, *(p+(i<<4)));
		// 	}
		// 	p[1] = 0;
		// 	_GLEXCOLL.SHM_buf =(volatile char*) p;
		// 	MPI_Barrier(Comm_intra);

		// }else{
		// 	//写入者
		// 	//其它进程只能打开共享缓存文件
		// 	MPI_Barrier(Comm_intra);
		// 	int fd = shm_open(host_name, O_RDWR, S_IRUSR|S_IWUSR);
		// 	ftruncate(fd,size);
		// 	char * pt = (char *)mmap(NULL,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
		// 	if (pt == MAP_FAILED)
		// 	{
		// 		printf("error map\n");
		// 		return ;
		// 	}
		// 	//SHM_BUF直接指向PT，计算过程中第I个块代表给I号发送消息
		// 	_GLEXCOLL.SHM_buf = pt;
		// }
		// //puts("init finished 1");
	}
	else if (strcmp(type, "RECURSIVE_DOUBLING") == 0)
	{
		// printf("%s\n",type);

		INTRA_ALLREDUCE_TYPE = 2;
		INIT_Shm_buf();
	}
	else if (strcmp(type, "L2_CACHE_AWARE") == 0)
	{
		// puts("check L2_CACHE_AWARE");
		INTRA_ALLREDUCE_TYPE = 3;
	}
	else if (strcmp(type, "OFF") == 0)
	{
		INTRA_ALLREDUCE_TYPE = -1;
	}

	MPI_Barrier(_GLEXCOLL.Comm_intra);
}

//获取一个K-nomial树的深度
int K_nominal_tree_stepN(int procn, int k)
{
	int stepn = 0;
	while (procn > 1)
	{
		//(procn/k)取上整
		// if(global_rank == 0)
		// 	printf("procn = %d\n",procn);
		procn = (procn + k - 1) / k;
		stepn++;
	}
	return stepn;
}

//获取自己在K-nomial树上的深度
int K_nominal_tree_my_stepN(int rank, int procn, int k)
{
	int stepn = 1;
	int stepN = K_nominal_tree_stepN(procn, k);
	for (stepn = 1; stepn < stepN; stepn++)
	{
		if (rank % k != 0)
			break;
		rank = rank / k;
	}
	return stepn;
}
//默认0为root
//第0步的parent为 rank - rank%k
//第1步的parent为 rank - rank%(k*2)
//第n-1不的parent为 rank - rank%(k ** (n))
int K_nominal_tree_parent(int rank, int k, int step)
{
	int K1 = k;
	for (int i = 0; i < step; i++)
		K1 *= k;
	return rank - rank % (K1);
}

void K_nominal_tree_child_vec(int rank, int procn, int step, int *Childvec, int *childn, int k)
{
	int K1 = 1;
	for (int i = 0; i < step; i++)
		K1 *= k;
	if (rank % (K1 * k) != 0)
	{
		*childn = -1;
		return;
	}

	int n = 0;
	for (n = 1; n < k; n++)
	{
		int tmp = rank + K1 * n;
		if (tmp >= procn)
			break;
		Childvec[n - 1] = tmp;
	}
	*childn = n - 1;
}

//大消息allreduce缓冲区的初始化

unsigned long allreduce_buffer_size_total = (1UL << 25);
unsigned long allreduce_buffer_size_single = (1UL << 16);
void *allreduce_intra_node_buf_cacheleaderBUF_header = 0;
void *large_allreduce_bufferP = 0;
void *large_allreduce_bufferROOT = 0;
int allreduce_CorePerCache = 4;
volatile void *allreduce_my_state_on_parent = 0;
volatile void *allreduce_cacheL_state_on_root = 0;
volatile void *allreduce_cacheL_flagsstart_on_root = 0;

volatile void *allreduce_shm_start;
volatile allreduce_Header *allreduce_shm_flags;
void *allreduce_shm_buffer_starts[64];
void Init_medium_allreduce_SH()
{
	unsigned long long allreduce_buffer_size_total = allreduce_buffer_size_single * intra_procn;
	if (intra_rank == 0)
	{
		char name[100];
		sprintf(name, "PJT-%s-Intra-large_msg_allreduce\n", host_name);
		int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
		if (fd < 0)
		{
			printf("shm_open failed!\n");
			return;
		}
		ftruncate64(fd, allreduce_buffer_size_total);
		allreduce_shm_start = mmap(NULL, allreduce_buffer_size_total, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		if (allreduce_intra_node_buf_cacheleaderBUF_header == MAP_FAILED)
		{
			puts("error map");
			exit(0);
		}
		memset((void *)allreduce_shm_start, 0, allreduce_buffer_size_total);
		allreduce_shm_flags = (volatile allreduce_Header *)allreduce_shm_start;
		allreduce_shm_buffer_starts[0] = (void *)allreduce_shm_start + sizeof(*allreduce_shm_flags);
		// *(allreduce_shm_flags->mlocks[10].lock0) = 996;
		MPI_Barrier(Comm_intra);
		// for (int i = 1; i < intra_procn;i++)
		// struct sync_header *
	}
	else
	{
		MPI_Barrier(Comm_intra);
		char name[100];
		sprintf(name, "PJT-%s-Intra-large_msg_allreduce\n", host_name);
		int fd = shm_open(name, O_RDWR, S_IRUSR | S_IWUSR);
		if (fd < 0)
		{
			printf("shm_open failed! 1\n");
			return;
		}
		allreduce_shm_start = mmap(NULL, allreduce_buffer_size_total, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		if (allreduce_intra_node_buf_cacheleaderBUF_header == MAP_FAILED)
		{
			puts("error map");
			exit(0);
		}
		allreduce_shm_flags = (volatile allreduce_Header *)allreduce_shm_start;
		// printf("        *(allreduce_shm_flags->mlocks[0].lock0) = %d\n",
		//        *(allreduce_shm_flags->mlocks[10].lock0));
		allreduce_shm_buffer_starts[0] = (void *)allreduce_shm_start + sizeof(*allreduce_shm_flags);
	}
	MPI_Barrier(Comm_intra);
	for (int i = 1; i < intra_procn; i++)
	{
		allreduce_shm_buffer_starts[i] = allreduce_shm_buffer_starts[i - 1] + allreduce_buffer_size_single;
	}
	MPI_Barrier(Comm_intra);
}

void init_allreduce_large_message_SH()
{
	if (intra_rank % allreduce_CorePerCache == 0)
	{
		char name[100];
		sprintf(name, "%s-%d-Intra-large_msg_allreduce\n", host_name, intra_rank);
		int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
		if (fd < 0)
		{
			printf("shm_open failed!\n");
			return;
		}
		ftruncate(fd, allreduce_buffer_size_total);
		allreduce_intra_node_buf_cacheleaderBUF_header = mmap(NULL, allreduce_buffer_size_total, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		if (allreduce_intra_node_buf_cacheleaderBUF_header == MAP_FAILED)
		{
			puts("error map");
			exit(0);
		}
		for (int i = 0; i < allreduce_buffer_size_total; i++)
		{
			// printf("i = %d\n", i);
			*(char *)(allreduce_intra_node_buf_cacheleaderBUF_header + i) = 0;
			// i % 26 + 'a';
		}
		struct sync_header *p = (struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header;
		for (int i = 0; i < allreduce_intra_node_locks_n; i++)
			*(volatile int *)(p->mlocks[i].lock0) = 0;
		large_allreduce_bufferP = allreduce_intra_node_buf_cacheleaderBUF_header + sizeof(allreduce_Header);
		// if (intra_rank == 0)
		//     printf("%p , %p\n", allreduce_intra_node_buf_cacheleaderBUF_header, large_allreduce_bufferP);
		MPI_Barrier(Comm_intra);

		if (intra_rank != 0)
		{
			//其它cache leader节点必须打开root节点的地址空间。
			sprintf(name, "%s-%d-Intra-large_msg_allreduce\n", host_name, 0);
			int fd = shm_open(name, O_RDWR, S_IRUSR | S_IWUSR);
			if (fd < 0)
			{
				printf("shm_open failed!\n");
				return;
			}
			int c = intra_rank / allreduce_CorePerCache;
			allreduce_cacheL_flagsstart_on_root = mmap(NULL, allreduce_buffer_size_total, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			allreduce_cacheL_state_on_root = ((struct sync_header *)allreduce_cacheL_flagsstart_on_root)->mlocks[c + allreduce_CorePerCache].lock0;
		}
	}
	else
	{
		MPI_Barrier(Comm_intra);
		int cache_leader = (intra_rank / allreduce_CorePerCache) * allreduce_CorePerCache;
		char name[100];
		sprintf(name, "%s-%d-Intra-large_msg_allreduce\n", host_name, cache_leader);
		int fd = shm_open(name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
		allreduce_intra_node_buf_cacheleaderBUF_header = mmap(NULL, allreduce_buffer_size_total, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		large_allreduce_bufferP = allreduce_intra_node_buf_cacheleaderBUF_header + sizeof(allreduce_Header);
		// if (intra_rank == 1)
		//     for (int i = 0; i < allreduce_buffer_size_total; i++)
		//     {
		//         printf("i = %c\n", *(char *)(allreduce_intra_node_buf_cacheleaderBUF_header + i));
		//     }
		struct sync_header *p = (struct sync_header *)allreduce_intra_node_buf_cacheleaderBUF_header;
		allreduce_my_state_on_parent = p->mlocks[intra_rank - cache_leader].lock0;
	}
	MPI_Barrier(Comm_intra);
}

double *allreduce_k_ary_bufvec[64];
char *allreduce_sendbuf, *allreduce_recvbuf;
int allreduce_send_recv_pair;
int allreduce_recursive_doubling_stepn;
int allreduce_buf_length;
MPI_Comm *vec_self_tuned_two_dimentional_allreduce_row_comm;
MPI_Comm *vec_self_tuned_two_dimentional_allreduce_colum_comm;
int vec_self_tuned_two_dimentional_allreduce_comm_NUM;
void GLEXCOLL_InitAllreduce()
{
	//创建共享通信缓冲区
	char *pathvar = getenv("INTRA_ALLREDUCE_TYPE");

	if (pathvar == 0)
	{
		puts("GLEXCOLL_InitAllreduce INTRA_ALLREDUCE_TYPE=0");
		exit(1);
	}
	if (ppn > 1)
		create_shemem(pathvar);
	// if (ppn % allreduce_CorePerCache == 0)
	// {
	// 	//创建节点内小中大消息共享内存
	// 	init_allreduce_large_message_SH();
	// }
	Init_medium_allreduce_SH();
	MPI_Barrier(MPI_COMM_WORLD);
	//  puts("1068");
	// MPI_Barrier(MPI_COMM_WORLD);

	//现在初始化规约树的每一层
	if (intra_rank == 0 && inter_procn > 1)
		_GLEXCOLL_Init_ReduceTree();
	if (intra_procn > 1)
	{
		//将TreeNumber的数值广播到节点内其它进程。
		MPI_Bcast(&TreeNumber, 1, MPI_INT, 0, Comm_intra);
	}
	// puts("1075");
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(Comm_inter);
	//  puts("1087");
	if (inter_procn > 1)
	{
		if (intra_rank == 0)
		{
			//初始化K-ary tree需要的缓冲区。
			for (int i = 0; i < 64; i++)
			{
				allreduce_k_ary_bufvec[i] = (double *)malloc(sizeof(double) * (1 << 13));
			}
			//接下来需要注册一对内存地址用于中小消息allreduce
			{
				allreduce_buf_length = (1 << 29);
				allreduce_sendbuf = (char *)malloc(allreduce_buf_length);
				allreduce_recvbuf = (char *)malloc(allreduce_buf_length);
				allreduce_send_recv_pair = GLECOLL_Register_SendRecv_Pair(allreduce_sendbuf, allreduce_buf_length,
																		  allreduce_recvbuf, allreduce_buf_length, Comm_inter);
			}

			int step = inter_procn;
			allreduce_recursive_doubling_stepn = 0;
			while (step > 1)
			{
				allreduce_recursive_doubling_stepn++;
				step = (step + 1) / 2;
			}

			// //初始化自适应二维划分通信子
			// if (inter_procn <= 2)
			// {
			// 	if (inter_rank == 0)
			// 		puts(" 节点数量太少，不能开启 自适应二维划分,不影响程序执行");
			// }
			// vec_self_tuned_two_dimentional_allreduce_row_comm = (MPI_Comm *)malloc(sizeof(*vec_self_tuned_two_dimentional_allreduce_row_comm) * inter_procn);
			// vec_self_tuned_two_dimentional_allreduce_colum_comm = (MPI_Comm *)malloc(sizeof(*vec_self_tuned_two_dimentional_allreduce_colum_comm) * inter_procn);
			// vec_self_tuned_two_dimentional_allreduce_comm_NUM = 0;
			// for (int dimX = 2; dimX <= inter_procn / 2; dimX++)
			// {
			// 	int x = inter_rank / dimX;
			// 	int y = inter_rank % dimX;
			// 	int color_row = x;
			// 	int color_colum = y;
			// 	MPI_Comm_split(Comm_inter, color_row, inter_rank, &(vec_self_tuned_two_dimentional_allreduce_row_comm[vec_self_tuned_two_dimentional_allreduce_comm_NUM]));
			// 	MPI_Comm_split(Comm_inter, color_colum, inter_rank, &(vec_self_tuned_two_dimentional_allreduce_colum_comm[vec_self_tuned_two_dimentional_allreduce_comm_NUM]));
			// 	vec_self_tuned_two_dimentional_allreduce_comm_NUM++;
			// }
		}
		else
		{
			allreduce_sendbuf = (char *)malloc(128);
			allreduce_recvbuf = (char *)malloc(128);
		}
	}
	else
	{
		allreduce_buf_length = (1 << 29);
		allreduce_sendbuf = (char *)malloc(allreduce_buf_length);
		allreduce_recvbuf = (char *)malloc(allreduce_buf_length);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(_GLEXCOLL.global_rank == 0)
	{
		puts("=======================GLEXCOLL_InitAllreduce初始化完成=======================");
		fflush(stdout);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// puts("1115");
}
