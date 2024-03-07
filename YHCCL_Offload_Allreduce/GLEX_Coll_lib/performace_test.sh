#!/bin/sh
echo "lets start"
for((noden=$MAX_NODEN;noden>=2;noden-=10))
do
    for((corePnuma=5;corePnuma<=5;corePnuma++))
    do
        # proc=16
        # corePnuma=0
        var=0
        let "var=1<<$corePnuma"
        let "proc=$noden*$GLEX_COLL_PPN"
        for((K=8;K>=2;K--))
        do
            for((calc=0;calc<$MAX_CALC;calc++))
            do
                echo "----------------------proc = $proc calc=$calc corePnuma =$var K = $K----------------------"
                cmd="yhrun -N $noden -n $proc ./build/test/$procName $calc 0 $corePnuma $K"
                echo $cmd
                for((loop=0;loop<10;loop++))
                do
                    yhrun -N $noden -n $proc -W 100 ./build/test/$procName $calc 0 $corePnuma $K
                    wait
                done
            done
        done
    done
done
#mpiexec  -n `expr $node * $GLEX_COLL_PPN` -ppn $GLEX_COLL_PPN ./build/test/main 0 0
#
