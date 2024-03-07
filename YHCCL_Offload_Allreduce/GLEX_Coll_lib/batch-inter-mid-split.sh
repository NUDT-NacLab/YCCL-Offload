#!/bin/sh


echo "-------------------------------------------------------------------------"
for((loop=0;loop<10;loop++))
do
#   节点间规约算法选择环境变量
    export  INTER_ALLREDUCE_TYPE=ORIGINAL
    yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName 0 0 5 2
    wait
    sleep 1
#   节点间规约算法选择环境变量
    export  INTER_ALLREDUCE_TYPE=PERFORMANCE_AWARENESS
    yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName 0 0 5 2
    wait
    sleep 1
    echo "-------------------------------------------------------------------------"
done 


# for((loop=0;loop<10;loop++))
# do
#    echo "-------------------------------loop =$loop-------------------------------------------------------------------------"
#    yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName 0 1
# done


# for((a=0;a<64;a++))
# do
#    for((b=0;b<64;b++))
#    do
#       if [ $a != $b ]
#       then
#          echo "-----------------------------------------------loop =$loop---------------------------------------------------------"
#          echo "yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName $a $b"
#          for((loop=0;loop<10;loop++))
#          do
#             yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName $a $b
#             wait
#          done
#       fi
#    done
# done

# for((a=0;a<16;a++))
# do
#    for((b=0;b<16;b++))
#    do
#       for((c=0;c<16;c++))
#       do
#          for((d=0;d<16;d++))
#          do
#             if [ $a != $b -a $c != $d ]
#             then
#                echo "-------------------------------start loop-------------------------------------------------------------------------"
#                echo "yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName $a $b $c $d"
#                for((loop=0;loop<5;loop++))
#                do
#                   yhrun -N $GLEX_COLL_NODEN -n $GLEX_COLL_PROCN ./build/test/$procName $a $b $c $d
#                   wait
#                done
#             fi
#          done
#       done
#    done
# done

