#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

struct Qdata{
    int start;
    int end;
};
typedef struct queue
{
    struct Qdata * pBase;//数组基址
    int front; //队首元素在数组中的下标
    int rear;  //队尾元素在数组中的下标
    int capacity;
}QUEUE;

void init_queue(QUEUE * pQ,int size);
int en_queue(QUEUE * pQ, int start,int end); //入队
int full_queue(QUEUE * pQ);
int out_queue(QUEUE * pQ,int *start,int *end); //出队
int empty_queue(QUEUE * pQ);

void free_queue(QUEUE * pQ);
#endif