#!/bin/bash
# export ALLTOALL_TYPE=MPI_LINEAR_EXCHANGE_COMPRESS
# for((leadn=4;leadn>=1;leadn-=1))
# do
#    for((ongoing=4;ongoing>=1;ongoing-=2))
#    do
#       echo "-------------------ongoing = 1<<$ongoing --leadn = $leadn----------------------------- "
#       yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName $leadn $ongoing
#       #yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_NODEN ./build/test/$procName $leadn $ongoing
#       wait
#       sleep 1
#    done 
# done

# export  ALLTOALL_TYPE=DIRECT
# yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_NODEN  ./build/test/$procName 4 32
# sleep 2
# export  ALLTOALL_TYPE=DIRECT_Kleader_NODE_AWARE_RDMA
# yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_NODEN  ./build/test/$procName 4 32

GLEX_COLL_NODEN_MAX=$GLEX_COLL_NODEN 
# if test -z "$Rack_ID";then
# 	echo "VARNAME is null string or VARNAME not exist"
# else
# 	echo $Rack_ID
#     RACK_FLAG=-M $Rack_ID
# fi


# for((noden=$GLEX_COLL_NODEN_MIN;noden<=$GLEX_COLL_NODEN_MAX;noden+=2))
# for((noden=$GLEX_COLL_NODEN_MAX;noden>=$GLEX_COLL_NODEN_MIN;noden/=2))
for((noden=$GLEX_COLL_NODEN_MIN;noden<=$GLEX_COLL_NODEN_MAX;noden*=2))
do

    export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $noden`
    # echo "------------node $noden procn $GLEX_COLL_PROCN msg_size=[0,1,2,...,20]*double ongoing=[1,2,4,8,14,32]------------"
    echo "---------------------node $noden procn $GLEX_COLL_PROCN------------------------------"
    # echo "------------------DIRECT----------------------------"
    # export  ALLTOALL_TYPE=DIRECT
    export NUMA_Aware_leaders=YES
    # export  ALLTOALL_TYPE=DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE
    yhrun -p $PARTITION -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 1 32
    # sleep 1.5
    # echo "------------------BRUCK----------------------------"
    # export  ALLTOALL_TYPE=BRUCK
    # yhrun -p $PARTITION -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 32
    # sleep 1.5
    # # echo "------------------MPI_LINEAR_EXCHANGE_COMPRESS----------------------------"
    # # export  ALLTOALL_TYPE=MPI_LINEAR_EXCHANGE_COMPRESS
    # # yhrun -p $PARTITION -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 32
    # # sleep 1.5
    # echo "------------------DIRECT_Kleader_NODE_AWARE_RDMA----------------------------"
    # export  ALLTOALL_TYPE=DIRECT_Kleader_NODE_AWARE_RDMA
    # yhrun -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 32
    # sleep 1
    # echo "------------------NUMA-Aware----------------------------"
    # export NUMA_Aware_leaders=YES
    # export  ALLTOALL_TYPE=DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE
    # echo " yhrun  -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 4"
    # yhrun  -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 4
    # sleep 1
    # echo "------------------Original----------------------------"
    # export NUMA_Aware_leaders=NO
    # export  ALLTOALL_TYPE=DIRECT_Kleader_NODE_AWARE_RDMA_PIPELINE
    # echo " yhrun  -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 4"
    # yhrun  -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 4
    # sleep 1
    # yhrun $RACK_FLAG -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName1 4 4
    # sleep 1
    # echo "------------------TPDS17_Cache_oblivious_intra_node----------------------------"
    # export  ALLTOALL_TYPE=TPDS17_Cache_oblivious_intra_node
    # yhrun -p $PARTITION  -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 4 32
    # sleep 1
done
# export noden=1
# export GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $noden`
# echo "------------node $noden1 msg_size=[0,1,2,...,20]*double ongoing=[1,2,4,8,14,32]------------"
# export  ALLTOALL_TYPE=DIRECT
# yhrun -p $PARTITION  -N $noden1 -n $GLEX_COLL_PROCN  ./build/test/$procName 4 32
# sleep 1.5
# export  ALLTOALL_TYPE=MPI_LINEAR_EXCHANGE_COMPRESS
# yhrun -p $PARTITION  -N $noden1 -n $GLEX_COLL_PROCN  ./build/test/$procName 4 32
# sleep 1.5