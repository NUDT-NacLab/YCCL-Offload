#!/bin/bash

#装载libnuma库，用于numa感知
dos2unix *.sh
export  GLEX_COLL_PPN=4
export  PARTITION=653mt2
#653
#test653
#test653
#th_mt1
#指定本机上所使用的python类型
export IPH_PYTHON=python
export  procName=allreduce-multipledata
#bcast-test
#alltoall-test-integer
#
#
#bcast-test
#
#alltoall-test
#
#alltoallv-test
#alltoallv-test
#
#
#alltoall-test
#allreduce-multipledata
cp /usr/local/glex/include/glex.h ./include
cp /usr/local/glex/lib/libglex.a ./lib
export  procName1=alltoall-test
export NUMA_Aware_leaders=NO
#alltoallv-test
#alltoall-test
#main-double
#alltoall-test
#intra-node-performance-test
#
#intra_node_contention_test
#intra-node-performance-test
#intra_node_contention_test
#intra_node_bigMSG
#
#
#intra_node_contention_test

export CUR_DIR=`pwd`
#节点间规约算法选择环境变量
export  INTER_ALLREDUCE_TYPE=ORIGINAL
#TJ_3_TOPOLOGY_AWARE
#PERFORMANCE_AWARENESS
#ORIGINAL
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
#BRUCKf
module add mpi-x
#module add mpi3-shared
export Rack_ID=
#cab0_4
#cab0_4 
#cab0_4
#FT cab0_4
#MT cab5_19
#MT cab50_59

#开启调试模式
# export PJT_GDB_DEBUG=1

cp /usr/local/glex/include/glex.h ./include/
#开启自动脚本取消
python ./job-auto-cancel.py $Rack_ID

# rm CMakeCache.txt
# cmake .
make clean

make glexcoll -j 16
# make glexcoll
make $procName
# make $procName1
rm slurm*
# exit


# exit 
nl=cn[3458-3866]
# cp $CUR_DIR/build/lib/* /vol8/home/nudtgcy/PJT-WorkDir/hpcg-3.1/glexcoll/
# cp $CUR_DIR/include/* /vol8/home/nudtgcy/PJT-WorkDir/hpcg-3.1/glexcoll/
rm alltoallv_performance.txt
rm alltoallv_all_statistics.txt
rm *.pdf

export EXCLUDE_NODE_LIST=cn1856,cn1860,cn1861,cn389,cn1929,cn1862
EXCLUDE_NODE_LIST_LENGTH=0
#24
EXCLUDE_NODE_LISTZNJ=
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
        export  GLEX_COLL_NODEN=1024
        export  GLEX_COLL_NODEN_MIN=512
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $GLEX_COLL_NODEN`
        Batch_NODEN=`expr $GLEX_COLL_NODEN + $EXCLUDE_NODE_LIST_LENGTH`
        echo "  python ../tools/TH3-NodeGenerator.py $Batch_NODEN $PARTITION $Rack_ID"
        # python ../tools/TH3-NodeGenerator.py $Batch_NODEN $PARTITION $Rack_ID
        #开启cab0_4 test653分区自动共享内存清除。需要root权限
        python3 ../tools/auto-memory-free.py nodefile.txt
        echo "yhbatch -M $Rack_ID -p $PARTITION  -N $Batch_NODEN --nodefile=nodefile.txt --time=0-0:15:20  ./batch-alltoall.sh"
        # exit
        # yhbatch -M $Rack_ID -p $PARTITION  -N $Batch_NODEN --nodefile=nodefile.txt --time=0-4:15:20  ./batch-alltoall.sh
        yhbatch -M $Rack_ID -p $PARTITION  --exclude=$EXCLUDE_NODE_LIST -N $GLEX_COLL_NODEN   --time=1-0:15:20  ./batch-alltoall.sh
        #  bash ./batch-alltoall.sh
        ;;
    8   )
        echo "checking specific nodelist"
        export  GLEX_COLL_NODEN=64
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $GLEX_COLL_NODEN`
        echo  echo \"  " `cat ./batch.sh | sed -n '2,2p'`" \" | bash
        echo "$EXCLUDE_NODE_LIST" > exclude_node_list.txt
        python GenNodeFile.py $GLEX_COLL_NODEN $PARTITION
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
        export  GLEX_COLL_NODEN_MAX=1024
        export  GLEX_COLL_NODEN_MIN=2
        wait
        yhbatch -p $PARTITION  --exclude=$EXCLUDE_NODE_LIST -N $GLEX_COLL_NODEN_MAX --time=2-6:33:20  ./batch-alltoallv.sh
         ;;
    11   )
        echo " allreduce small mid big message Testing"
        export  GLEX_COLL_NODEN_MAX=200
        export  GLEX_COLL_NODEN_MIN=200
        wait
        # cn[20992-21247]
        # cn[23648-23655]
        # cn23648,cn23649,cn23651,cn23654
        #  
        # -w cn[20984-21239]

        mcmd="yhbatch -p $PARTITION  -x `cat ./exclude-th3.txt` -N $GLEX_COLL_NODEN_MAX --time=2-4:20:20  ./batch-allreduce.sh"
        `echo $mcmd`
        # yhbatch -p $PARTITION -w cn23648,cn23649,cn23651,cn23654 -N $GLEX_COLL_NODEN_MAX --time=2-4:20:20  ./batch-allreduce.sh
         ;;
    12   )
        echo " Broadcast Testing"
        export  GLEX_COLL_NODEN_MAX=8
        export  GLEX_COLL_NODEN_MIN=8
        Batch_NODEN=`expr $GLEX_COLL_NODEN_MAX + $EXCLUDE_NODE_LIST_LENGTH`
        python GenNodeFile.py $Batch_NODEN  $PARTITION  
        # yhrun -p $PARTITION --nodefile=nodefile.txt -N $GLEX_COLL_NODEN_MAX -n $GLEX_COLL_NODEN_MAX pwd
        wait
        # yhbatch -p $PARTITION -N $Batch_NODEN  -F hostfile.txt --time=0-2:13:20  ./batch-broadcast.sh
        yhbatch -p $PARTITION  --exclude=$EXCLUDE_NODE_LIST -N $Batch_NODEN --time=0-1:3:20  ./batch-broadcast.sh
        # bash ./batch-broadcast.sh > out.txt
         ;;

esac
for((i=1;i<=1000;i++))
do
    echo "---------------------------------------------------------------------------------------"
    cat slurm* 
    sleep 2
    echo "---------------------------------------------------------------------------------------"
done