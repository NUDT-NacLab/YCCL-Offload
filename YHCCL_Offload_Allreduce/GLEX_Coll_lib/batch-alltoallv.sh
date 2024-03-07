#!/bin/bash
# for((leadn=6;leadn>=1;leadn-=1))
# do
#    for((ongoing=4;ongoing>=1;ongoing-=2))
#    do
#       echo "-------------------ongoing = 1<<$ongoing --leadn = $leadn----------------------------- "
#       #yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName $leadn $ongoing
#       yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_NODEN ./build/test/$procName $leadn $ongoing
#       wait
#       sleep 1
#    done 
# done
##for((noden=$GLEX_COLL_NODEN_MIN;noden<=$GLEX_COLL_NODEN_MAX;noden+=50))
for((noden=$GLEX_COLL_NODEN_MIN;noden<=$GLEX_COLL_NODEN_MAX;noden*=2))
do
    for((ongoing=1;ongoing>=1;ongoing-=1))
    do
        export  GLEX_COLL_PROCN=`expr $GLEX_COLL_PPN \* $noden`
        # for((loop=0;loop<1;loop++))
        for((sz=10;sz<=25;sz+=1))
        do
            echo "--------noden=$noden, procn=$GLEX_COLL_PROCN, ongoing=$ongoing sz=1<<$sz---------------------------------"
            export AlltoallvAlgorithmTYPE=LINEAR_SHIFT
            yhrun -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName $ongoing $sz
            sleep 1
            wait
            echo ""
            # mv ./alltoallv_all_statistics.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-step-performance.py $CUR_DIR/data_visit/  "$AlltoallvAlgorithmTYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait
            # mv ./alltoallv_performance.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-min-avg-max-steptime.py $CUR_DIR/data_visit/  "stepinfo-$AlltoallvAlgorithmTYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait


            # # export AlltoallvAlgorithmTYPE=XOR_EXCHANGE
            # # yhrun -N $noden -n $noden  ./build/test/$procName $ongoing 1
            # # sleep 1
            # # wait

            # export AlltoallvAlgorithmTYPE=OPT_BANDWIDTH
            # yhrun -N $noden -n $noden  ./build/test/$procName 1 1
            # sleep 1
            # wait
            # echo""
            # yhrun -N $noden -n $noden  ./build/test/$procName 2 1
            # sleep 1
            # wait
            # echo""
            
            # export AlltoallvAlgorithmTYPE=LINEAR_SHIFT_REORDER
            # export ORDER_TYPE=reverse_order
            # yhrun -N $noden -n $noden  ./build/test/$procName $ongoing 1
            # sleep 1
            # mv ./alltoallv_all_statistics.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-step-performance.py $CUR_DIR/data_visit/  "$AlltoallvAlgorithmTYPE-$reverse_order-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait
            # mv ./alltoallv_performance.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-min-avg-max-steptime.py $CUR_DIR/data_visit/  "stepinfo-$AlltoallvAlgorithmTYPE--$reverse_order-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait
            
            
            # export AlltoallvAlgorithmTYPE=LINEAR_SHIFT_REORDER
            # export ORDER_TYPE=random_order
            # yhrun -N $noden -n $noden  ./build/test/$procName $ongoing 1
            # sleep 1
            # # mv ./alltoallv_all_statistics.txt ./data_visit
            # # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-step-performance.py $CUR_DIR/data_visit/  "$AlltoallvAlgorithmTYPE-$ORDER_TYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # # wait
            # mv ./alltoallv_performance.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-min-avg-max-steptime.py $CUR_DIR/data_visit/  "stepinfo-$AlltoallvAlgorithmTYPE--$ORDER_TYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait

            # export AlltoallvAlgorithmTYPE=LINEAR_SHIFT_REORDER
            # export ORDER_TYPE=diavation_order
            # yhrun -N $noden -n $noden  ./build/test/$procName $ongoing 1
            # sleep 1
            # # mv ./alltoallv_all_statistics.txt ./data_visit
            # # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-step-performance.py $CUR_DIR/data_visit/  "$AlltoallvAlgorithmTYPE-$ORDER_TYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # # wait
            # mv ./alltoallv_performance.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-min-avg-max-steptime.py $CUR_DIR/data_visit/  "stepinfo-$AlltoallvAlgorithmTYPE--$ORDER_TYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait

            # export AlltoallvAlgorithmTYPE=LINEAR_SHIFT_REORDER
            # export ORDER_TYPE=ring_distance_order
            # yhrun -N $noden -n $noden  ./build/test/$procName $ongoing 1
            # sleep 1
            # mv ./alltoallv_all_statistics.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-step-performance.py $CUR_DIR/data_visit/  "$AlltoallvAlgorithmTYPE-$ORDER_TYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait
            # mv ./alltoallv_performance.txt ./data_visit
            # $IPH_PYTHON $CUR_DIR/data_visit/alltoallv-min-avg-max-steptime.py $CUR_DIR/data_visit/  "stepinfo-$AlltoallvAlgorithmTYPE--$ORDER_TYPE-noden=$noden-procn=$noden-ongoing=$ongoing"
            # wait
            
            # 使用每节点单进程
            # export AlltoallvAlgorithmTYPE=RDMA_LINEAR_SHIFT_REORDER
            # export ORDER_TYPE=random_order
            # # export AlltoallvAlgorithmTYPE=LINEAR_SHIFT
            # yhrun -p $PARTITION -N $noden -n $noden  ./build/test/$procName1 $ongoing 1
            # sleep 2

            # #没节点32进程
            # export AlltoallvAlgorithmTYPE=LINEAR_SHIFT_COMPRESS
            # # echo "yhrun -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName $ongoing 1"
            # yhrun -N $noden -n $GLEX_COLL_PROCN  ./build/test/$procName $ongoing 1
            # # yhrun -p all -N $noden -n 4  pwd
            # sleep 1
            # wait

            wait
            echo "------------------------------------"
        done
    done
done

