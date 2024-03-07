#!/bin/bash
cp /usr/local/glex/include/glex.h ./include
cp /usr/local/glex/lib/libglex.a ./lib
#装载libnuma库，用于numa感知
dos2unix *.sh
export  GLEX_COLL_PPN=4
export  PARTITION=bigdata
#
#th_mt1
#指定本机上所使用的python类型
export IPH_PYTHON=python
#
#all
#th_mt1
#all
# module add gcc/7.5.0
# module load cmake/3.20.2-gcc-7.5.0

#
export  procName=allreduce-multipledata
export NUMA_Aware_leaders=NO
#alltoall-test-integer
#bcast-test
#alltoall-test
#
#alltoallv-test
#alltoallv-test
#
#
#alltoall-test
#allreduce-multipledata
#
#
export  procName1=alltoall-test
#alltoallv-test
#alltoall-test
#main-double
#alltoall-test
#intra-node-performance-test
#
#intra_node_contention_test
#intra-node-performance-test
#intra_node_contention_test
#
#
#intra_node_bigMSG
#
#
#intra_node_contention_test

export CUR_DIR=`pwd`
#节点间规约算法选择环境变量
export  INTER_ALLREDUCE_TYPE=ORIGINAL
#PERFORMANCE_AWARENESS
#ORIGINAL
#
#
#PERFORMANCE_AWARENESS
#节点内规约算法选择环境变量
export  INTRA_ALLREDUCE_TYPE=L2_CACHE_AWARE

#REDUCE_BCAST
#L2_CACHE_AWARE
#REDUCE_BCAST
#REDUCE_BCAST_REORDER
#
#L2_CACHE_AWARE
#RECURSIVE_DOUBLING
#2_NOMINAL
#
#REDUCE_BCAST
#
#选择alltoall算法
export  ALLTOALL_TYPE=DIRECT_Kleader_NODE_AWARE_RDMA
#BRUCK_RDMA
#
#DIRECT_Kleader_NODE_AWARE_RDMA
#DIRECT_Kleader_NODE_AWARE
#DIRECT_NODE_AWARE_alltoall
#DIRECT_NODE_AWARE_alltoall
#DIRECT_NODE_AWARE_RDMA
#DIRECT_NODE_AWARE
#
#BRUCK
#module add MPI/mpich3.2.1/gcc8.2.0
#module add MPI/mpich3.2.1/gcc8.2.0-noglex
export LD_LIBRARY_PATH=/BIGDATA1/nudt_jliu_1/pjt/yh-coll/GLEX_Coll_lib:$LD_LIBRARY_PATH
cp /usr/local/glex/include/glex.h ./include/
module add cmake/3.20.2


#开启调试模式
# export PJT_GDB_DEBUG=1

#开启自动脚本取消
python ./job-auto-cancel.py
# rm CMakeCache.txt
# cmake .
make clean
make glexcoll -j 8
# make glexcoll
make $procName
# make $procName1
rm slurm*

nl=cn[3458-3866]
rm alltoallv_performance.txt
rm alltoallv_all_statistics.txt
rm *.pdf

EXCLUDE_NODE_LIST=
EXCLUDE_NODE_LISTZNJ=
EXIT 0
MODE=11
case $MODE in
    1   )
        echo "checking  correctiveness"
        export  GLEX_COLL_NODEN=512
        export  GLEX_COLL_NODENMin=512
        export  GLEX_COLL_PROCN=`expr 32 \* $GLEX_COLL_NODEN`
        echo  echo \"  " `cat ./batch.sh | sed -n '2,2p'`" \" | bash
        yhbatch -p $PARTITION --exclude=$EXCLUDE_NODE_LIST  -N $GLEX_COLL_NODEN --time=0-0:14:20  ./batch.sh
        ;;
    2   )
        echo "performance testing"
        export MAX_NODEN=30
        export MAX_PROC=960
        yhbatch -p $PARTITION --exclude=$EXCLUDE_NODE_LIST -N $MAX_NODEN --time=3-3:1:20  ./performace_test.sh
        ;;
    3   )
        echo "nonblocking performance testing"
        export MAX_NODEN=512
        export MAX_PROC=16384
        export MAX_CALC=1
        yhbatch -p $PARTITION --exclude=$EXCLUDE_NODE_LIST -N $MAX_NODEN -n $MAX_PROC --time=4-3:1:20  ./performace_test.sh
        ;;
    4   )
        echo "intra algorithm testing"
        export  GLEX_COLL_PROCN=32
        export  GLEX_COLL_NODEN=1
        yhbatch -p $PARTITION --exclude=$EXCLUDE_NODE_LIST -N $GLEX_COLL_NODEN --time=0-0:5:20  ./batch_intra.sh
        ;;
    5   )

        echo "inter algorithm testing"
        export  GLEX_COLL_PROCN=16384
        export  GLEX_COLL_NODEN=512
        yhbatch -p $PARTITION --exclude=$EXCLUDE_NODE_LIST -N $GLEX_COLL_NODEN --time=0-0:15:20  ./batch-inter-mid-split.sh
        ;;
        
    6   )
        echo "Start All Testing"
        export  GLEX_COLL_NODEN_END=20
        export  GLEX_COLL_NODEN_START=1
        export  MAX_CALC=10
        yhbatch -p $PARTITION --exclude=$EXCLUDE_NODE_LIST -N $GLEX_COLL_NODEN_END --time=0-9:15:20  ./batch-total.sh
        ;;

    7   )
        echo "alltoall Testing"
        export  GLEX_COLL_NODEN=48
        export  GLEX_COLL_NODEN_MIN=3
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $GLEX_COLL_NODEN`
        echo  echo \"  " `cat ./batch.sh | sed -n '2,2p'`" \" | bash
        yhbatch -p $PARTITION   -N $GLEX_COLL_NODEN --time=0-1:6:20  ./batch-alltoall.sh
        ;;
    8   )
        echo "checking specific nodelist"
        export  GLEX_COLL_NODEN=64
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $GLEX_COLL_NODEN`
        echo  echo \"  " `cat ./batch.sh | sed -n '2,2p'`" \" | bash
        python GenNodeFile.py $GLEX_COLL_NODEN 2
        wait
        yhbatch -p $PARTITION  --exclude=$EXCLUDE_NODE_LIST --nodefile=nodefile.txt -N $GLEX_COLL_NODEN --time=0-0:1:20  ./batch.sh
        ;;
    9   )
        echo "ZNJ alltoall Testing"
        export  GLEX_COLL_NODEN=3
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $GLEX_COLL_NODEN`
        wait
        yhbatch -p $PARTITION  --exclude=$EXCLUDE_NODE_LISTZNJ -N $GLEX_COLL_NODEN --time=0-0:1:20  ./batch-alltoall.sh
        ;;
    10   )
        echo " alltoallv Testing"
        export  GLEX_COLL_NODEN_MAX=128
        export  GLEX_COLL_NODEN_MIN=2
        wait
        yhbatch -p $PARTITION  --exclude=$EXCLUDE_NODE_LIST -N $GLEX_COLL_NODEN_MAX --time=1-6:33:20  ./batch-alltoallv.sh
         ;;
    11   )
        echo " allreduce small mid big message Testing"
        export  GLEX_COLL_NODEN_MAX=20
        export  GLEX_COLL_NODEN_MIN=20
        wait
        yhbatch -p $PARTITION -x `cat ./node-screening/excludeth2` -N $GLEX_COLL_NODEN_MAX --time=0-0:20:20  ./batch-allreduce.sh
         ;;
    12   )
        echo " Broadcast Testing"
        export  GLEX_COLL_NODEN_MAX=32
        export  GLEX_COLL_NODEN_MIN=4
        wait
        yhbatch -p $PARTITION  --exclude=$EXCLUDE_NODE_LIST -N $GLEX_COLL_NODEN_MAX --time=0-0:13:20  ./batch-broadcast.sh
         ;;

esac
for((i=1;i<=1000;i++))
do
    echo "---------------------------------------------------------------------------------------"
    cat slurm* 
    sleep 2
    echo "---------------------------------------------------------------------------------------"
done
