#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double tmp = 1.0;
    double target = 0.0;
    MPI_Request req;
    MPI_Status status;
    int loopN = 10000000;
    for (int i = 0; i < loopN; i++)
    {
        MPI_Iallreduce(&tmp, &target, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req);
        MPI_Wait(&req, &status);
        if (i % 10000 == 0)
            if (rank == 0)
                printf("i = %d\n", i);
    }

    MPI_Finalize();
    return 0;
}
