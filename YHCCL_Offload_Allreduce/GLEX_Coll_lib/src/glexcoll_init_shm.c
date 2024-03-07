
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <pthread.h>
#include <mpi.h>
#include "glexcoll.h"

static int buffersize = 1 << 13; //每个缓冲区8kb
volatile int *shared_release_flags[33];
volatile int *shared_gather_flags[33];
volatile char *broad_cast_buffer;
volatile char *Reduce_buffers[33];

int Get_shared_buffer_size()
{
    return 64 * 32 + 64 * 32 + buffersize + buffersize * 32;
}
void Intel_Init_SHMBUF()
{
    int size = Get_shared_buffer_size();
    sprintf(_GLEXCOLL.name, "%s-intel\0", host_name);
    volatile char *p;
    if (intra_rank == 0)
    {
        int fd = shm_open(_GLEXCOLL.name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
        ftruncate(fd, size);
        p = (volatile char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (p == MAP_FAILED)
        {
            printf("error map\n");
            return;
        }
        for (int i = 0; i < size; i++)
        {
            p[i] = 0;
        }
        MPI_Barrier(Comm_intra);
    }
    else
    {
        MPI_Barrier(Comm_intra);
        int fd = shm_open(_GLEXCOLL.name, O_RDWR, S_IRUSR | S_IWUSR);
        ftruncate(fd, size);
        p = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    }
    for (int i = 0; i < intra_procn; ++i)
    {

        //printf("intra_procn = %d\n",intra_procn);
        shared_release_flags[i] = (volatile int *)p;
        p = p + 64;
    }
    for (int i = 0; i < intra_procn; ++i)
    {
        shared_gather_flags[i] = (volatile int *)p;
        p += 64;
    }
    broad_cast_buffer = p;
    p += buffersize;
    for (int i = 0; i < intra_procn; ++i)
    {
        Reduce_buffers[i] = p;
        p += buffersize;
    }
}
void Intel_release(int bcastflag, void *sendbuf, void *recvbuf, int count)
{
    static int pvt_release_state = 0;
    if (bcastflag == 1 && intra_rank == 0)
    {
        for (int c = 0; c < count; c++)
        {
            ((double *)broad_cast_buffer)[c] = ((double *)sendbuf)[c];
        }
    }
    __sync_synchronize();
    pvt_release_state++;
    if (intra_rank == 0)
    {
        //rooot
        *(shared_release_flags[0]) = pvt_release_state;
        if (bcastflag)
            for (int c = 0; c < count; c++)
            {
                ((double *)recvbuf)[c] = ((double *)sendbuf)[c];
            }
    }
    else
    {
        //non-root
        while (*(shared_release_flags[0]) != pvt_release_state)
            ;
        if (bcastflag)
        {
            __sync_synchronize();
            for (int c = 0; c < count; c++)
            {
                ((double *)recvbuf)[c] = ((double *)broad_cast_buffer)[c];
            }
        }
    }
}
void Intel_gather(int reduceflag, void *sendbuf, void *recvbuf, int count)
{
    static int pvt_gather_state = 0;
    if (reduceflag)
    {
        if (intra_rank == 0)
        {
            for (int i = 0; i < count; i++)
            {
                ((double *)recvbuf)[i] = ((double *)sendbuf)[i];
            }
        }
        else
        {
            for (int i = 0; i < count; i++)
            {
                ((double *)Reduce_buffers[intra_rank])[i] = ((double *)sendbuf)[i];
            }
        }
    }
    pvt_gather_state++;
    if (intra_rank == 0)
    {
        for (int i = 1; i < intra_procn; ++i)
        {
            while (*(shared_gather_flags[i]) != pvt_gather_state)
                ;

            if (reduceflag)
            {
                __sync_synchronize();
                for (int c = 0; c < count; c++)
                    ((double *)recvbuf)[c] += ((double *)Reduce_buffers[i])[c];
            }
        }
        // if(reduceflag){
        //     printf("recvbuf[0] = %f\n",((double*)recvbuf)[0]);
        // }
    }
    __sync_synchronize();
    *(shared_gather_flags[intra_rank]) = pvt_gather_state;
}
// void Init_SHMBUF()
// {
//     int size = intra_procn * bufN_between_ab * bufsize;
//     char name[100];
// 	sprintf(name,"%s-%d\0",host_name,intra_rank);
// 	int fd = shm_open(_GLEXCOLL.name,O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
// 	ftruncate(fd,size);
//     volatile char * p = (char *)mmap(NULL,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
// 	if (p == MAP_FAILED)
// 	{
// 		printf("error map\n");
// 		return ;
// 	}
// 	memset((void *)p,'S',size);
// 	MPI_Barrier(Comm_intra);
// 	//接下来记录下每对进程之间共享内存的区域
// 	for(int i = 0;i<intra_procn;++i)
// 	{
// 		sprintf(name,"%s-%d\0",host_name,i);
// 		fd = shm_open(_GLEXCOLL.name, O_RDWR, S_IRUSR | S_IWUSR);
// 		ftruncate(fd,size);
// 		SHM_bufss[i] = mmap(NULL,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
// 	}
// 	MPI_Barrier(Comm_intra);
// }
// void Init_SHEMFLAG()
// {
//     int size = sizeof(struct SHM_flag_between_source_dest)*intra_procn;
//     char name[100];
// 	sprintf(name,"%s-%d-flag\0",host_name,intra_rank);
// 	int fd = shm_open(_GLEXCOLL.name,O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
// 	ftruncate(fd,size);
//     volatile struct SHM_flag_between_source_dest *p = (char *)mmap(NULL,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
// 	if (p == MAP_FAILED)
// 	{
// 		printf("error map\n");
// 		return ;
// 	}
//     for(int i = 0;i<intra_procn;++i)
//         for(int j = 0;j<bufN_between_ab;j++)
//             {
//                 p[i].flags[j].flag = 'I';
//                 p[i].flags[j].size = 0;
//             }
// 	MPI_Barrier(Comm_intra);
// 	//接下来记录下每对进程之间共享内存的区域
// 	for(int i = 0;i<intra_procn;++i)
// 	{
// 		sprintf(name,"%s-%d-flag\0",host_name,i);
// 		fd = shm_open(_GLEXCOLL.name, O_RDWR, S_IRUSR | S_IWUSR);
// 		ftruncate(fd,size);
// 		SHM_flags[i] = mmap(NULL,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
//         for(int j = 0;j<intra_procn;++j)
//             Using_buf_id[i][j] = 0;
// 	}

// 	MPI_Barrier(Comm_intra);
// }

// volatile char * Get_buf_to(int dest)
// {
//     return SHM_bufs[dest] + intra_rank * bufN_between_ab * bufsize + Using_buf_id[intra_rank][dest]*bufsize;
// }
// volatile char * Get_buf_from(int source)
// {
//     return SHM_bufs[intra_rank] + source * bufN_between_ab * bufsize + Using_buf_id[source][intra_rank]*bufsize;
// }
// void Set_flag_to_send(int dest,int size)
// {
//     SHM_flags[dest].flags[Using_buf_id[intra_rank][dest]].flag = 'S';
//     SHM_flags[dest].flags[Using_buf_id[intra_rank][dest]].size = size;
// }
// void Set_flag_to_recv(int source,int size)
// {
//     SHM_flags[source].flags[Using_buf_id[intra_rank][dest]].flag = 'S';
//     SHM_flags[source].flags[Using_buf_id[intra_rank][dest]].size = size;
// }
// void Send_Msg_to(int dest,char * sendbuf,int size)
// {

// }
// void Get_Msg_from(int source,char *recvbuf,int size)
// {

// }