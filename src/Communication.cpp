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