#!/bin/sh


for((noden=$GLEX_COLL_NODEN_END;noden>=$GLEX_COLL_NODEN_START;noden-=10))
do
    corePnuma=5
    # corePnuma=0
    var=0
    let "var=1<<$corePnuma"
    let "proc=$noden*$GLEX_COLL_PPN"
        echo "----------------------loop $loop noden $noden-------------------------"

        for((calc=0;calc<1;calc++))
        do
            echo "---------------------------MPI result--calc $calc------------------------------"
            for((K=8;K>=2;K--))
            do
                echo "----------------------proc = $proc calc=$calc corePnuma =$var K = $K----------------------"
                cmd="yhrun -N $noden -n $proc ./build/test/$procName $calc 0 $corePnuma $K"
                echo $cmd
                
                for((loop=0;loop<20;loop++))
                do
                    #测试两种规约树。
                    #第一个列为原始规约树
                    export  INTER_ALLREDUCE_TYPE=ORIGINAL
                    #ORIGINAL
                    #PERFORMANCE_AWARENESS
                    #数据测GLEX
                    yhrun -N $noden -n $proc -W 100 ./build/test/$procName $calc 0 $corePnuma $K 3
                    wait
                    # #第二列为距离平衡规约树
                    # export  INTER_ALLREDUCE_TYPE=PERFORMANCE_AWARENESS
                    # #ORIGINAL
                    # #PERFORMANCE_AWARENESS
                    # yhrun -N $noden -n $proc -W 100 ./build/test/$procName $calc 0 $corePnuma $K 1
                    # wait
                    echo ""
                done
            done
        done

done