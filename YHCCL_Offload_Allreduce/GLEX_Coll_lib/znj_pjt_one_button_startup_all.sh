#!/bin/bash
dos2unix *.sh
#装载libnuma库，用于numa感知
export  GLEX_COLL_PPN=32
export  PARTITION=all
export  procName=alltoall-test-integer
#alltoallv-test
export  procName1=alltoallv-rdma
#alltoallv-test
export CUR_DIR=`pwd`
#节点间规约算法选择环境变量
export  INTER_ALLREDUCE_TYPE=PERFORMANCE_AWARENESS
#指定本机上所使用的python类型
export IPH_PYTHON=python3
#ORIGINAL
#
#
#PERFORMANCE_AWARENESS
#节点内规约算法选择环境变量
export  INTRA_ALLREDUCE_TYPE=L2_CACHE_AWARE

#选择alltoall算法
export  ALLTOALL_TYPE=DIRECT_Kleader_NODE_AWARE_RDMA
#module add MPI/mpich-3.2.1-yh202005
#cp /usr/local/glex/include/glex.h ./include/
# dos2unix *.sh
# make clean
# rm CMakeCache.txt
# cmake .
# make clean
# make -j 4
make
rm slurm*
rm alltoallv_performance.txt
rm alltoallv_all_statistics.txt
rm *.pdf
python ./job-auto-cancel.py
nl=cn[3458-3866]
EXCLUDE_NODE_LIST=cn4136,cn4166
EXCLUDE_NODE_LISTZNJ=fn128
MODE=7
case $MODE in
    1   )
        echo "checking  correctiveness"
        export  GLEX_COLL_NODEN=40
        export  GLEX_COLL_NODENMin=40
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
        export  GLEX_COLL_NODEN=1
        module add tcl/tcl-8.5.19-gcc4.9.3
        #tcl/tcl-8.5.19-gcc8.2.0
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $GLEX_COLL_NODEN`
        # export  GLEX_COLL_PROCN=2048  
        echo  echo \"  " `cat ./batch.sh | sed -n '2,2p'`" \" | bash
        yhbatch -p $PARTITION  -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN --time=0-1:15:20  ./batch-alltoall.sh
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
        # module add GCC/GCC-8.3.0
        module add tcl/tcl-8.5.19-gcc8.2.0
        #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gv3/software/GCC/gcc-7.3.0/gcc730/lib64/
        export  GLEX_COLL_NODEN_MAX=512
        export  GLEX_COLL_NODEN_MIN=2
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $GLEX_COLL_NODEN_MAX`
        wait
        yhbatch -p $PARTITION  --exclude=$EXCLUDE_NODE_LISTZNJ -N $GLEX_COLL_NODEN_MAX -n $GLEX_COLL_PROCN --time=0-14:2:20  ./batch-alltoallv.sh
         ;;

esac
for((i=1;i<=1000;i++))
do
    echo "---------------------------------------------------------------------------------------"
    cat slurm* 
    sleep 2
    echo "---------------------------------------------------------------------------------------"
done