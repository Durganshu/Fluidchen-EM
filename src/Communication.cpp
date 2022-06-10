#include "Communication.hpp"

static void communicate(Fields &field, Domain &domain) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (domain.neighbours[0] != -1) {
        // Communicate to Left
    }

    if (domain.neighbours[1] != -1) {
        // Communicate to Right
    }

    if (domain.neighbours[2] != -1) {
        // Communicate to Top
    }

    if (domain.neighbours[3] != -1) {
        // Communicate to Bottom
    }
}

double reduce_min(double dt)
{
    int rank;int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double final_dt;

    MPI_Reduce(&dt,&final_dt,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Bcast(&final_dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    return final_dt;

}