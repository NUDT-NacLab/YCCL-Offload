#!/bin/bash
echo "lets start"
export  INTER_ALLREDUCE_TYPE=ORIGINAL
# sleep 3
# for((noden=$GLEX_COLL_NODEN_MAX;noden>=$GLEX_COLL_NODEN_MIN;noden/=2))
for((noden=$GLEX_COLL_NODEN_MIN;noden<=$GLEX_COLL_NODEN_MAX;noden*=2))
# for((noden=$GLEX_COLL_NODEN_MAX;noden>=$GLEX_COLL_NODEN_MIN;noden-=10))
do
    # noden=`expr 1 + $noden`
    export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $noden`
    # for((k=3;k>=3;k-=2))

export  INTER_ALLREDUCE_TYPE=TJ_3_TOPOLOGY_AWARE
            echo "-------------------------规约算法为:$INTER_ALLREDUCE_TYPE-----noden=$noden  K=$k-----------------------------------"
            yhrun --mpi=pmix -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 9 4
            wait
            sleep 4
    for((k=2;k<=8;k+=2))
    do
        for((loopn=0;loopn<1;loopn++))
        do 
export  INTER_ALLREDUCE_TYPE=ORIGINAL
            echo "-------------------------规约算法为:$INTER_ALLREDUCE_TYPE-----noden=$noden  K=$k-----------------------------------"
            yhrun --mpi=pmix -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName 9 $k
            wait
            sleep 4
        done
    done
done
exit 