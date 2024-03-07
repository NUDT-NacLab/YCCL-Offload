#!/bin/bash
yhrun -p $PARTITION -N 2 -n 2 pwd
# for((noden=$GLEX_COLL_NODEN_MIN;noden<=$GLEX_COLL_NODEN_MAX;noden*=2))
# for((noden=$GLEX_COLL_NODEN_MAX;noden>=$GLEX_COLL_NODEN_MIN;noden-=10))
for((noden=$GLEX_COLL_NODEN_MIN;noden<=$GLEX_COLL_NODEN_MAX;noden*=2))
do
    for((k=2;k<=2;k+=2))
    do
        for((loopn=0;loopn<1;loopn++))
        do 
            echo "------------------------------------noden=$noden  K=$k-------------------------------------------"
            yhrun -p $PARTITION  -N $noden -n $noden ./build/test/$procName 13 $k
            wait
            sleep 1.5
        done
    done
done