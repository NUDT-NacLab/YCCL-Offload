#include <stdlib.h>
#include "BasicDataStructure.h"

void init_queue(QUEUE *pQ, int size);
int en_queue(QUEUE *pQ, int start, int end); //入队
int full_queue(QUEUE *pQ);
int out_queue(QUEUE *pQ, int *start, int *end); //出队
int empty_queue(QUEUE *pQ);

void free_queue(QUEUE *pQ);

void init_queue(QUEUE *pQ, int size)
{
    pQ->pBase = (struct Qdata *)malloc(sizeof(struct Qdata) * size);
    pQ->front = 0;
    pQ->rear = 0;
    pQ->capacity = size;
}

int full_queue(QUEUE *pQ)
{
    if ((pQ->rear + 1) % (pQ->capacity) == pQ->front)
        return 1;
    else
        return 0;
}

int en_queue(QUEUE *pQ, int start, int end)
{
    if (full_queue(pQ))
        return 0;
    else
    {
        pQ->pBase[pQ->rear].start = start;
        pQ->pBase[pQ->rear].end = end;
        pQ->rear = (pQ->rear + 1) % (pQ->capacity);

        return 1;
    }
}

int empty_queue(QUEUE *pQ)
{
    if (pQ->front == pQ->rear)
        return 1;
    else
        return 0;
}

int out_queue(QUEUE *pQ, int *start, int *end)
{
    if (empty_queue(pQ)) //判断队列是否为空
        return 0;
    else
    {
        *start = pQ->pBase[pQ->front].start;
        *end = pQ->pBase[pQ->front].end;
        pQ->front = (pQ->front + 1) % (pQ->capacity);

        return 1;
    }
}
void free_queue(QUEUE *pQ)
{
    free(pQ->pBase);
}