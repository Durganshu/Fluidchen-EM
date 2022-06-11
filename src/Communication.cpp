#include "Communication.hpp"
#include <vector>

static void communicate(Matrix<double> &data, Domain &domain) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    // Communicate to Left
    if (domain.neighbours[0] != -1) {
        std::vector<double> send = data.get_col(1);
        std::vector<double> recv;
        MPI_Sendrecv(&send, domain.size_y + 2, MPI_DOUBLE, domain.neighbours[0], 1, &recv, domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[0], 2, MPI_COMM_WORLD, &status);
        data.set_col(recv, 0);
    }

    // Communicate to Right
    if (domain.neighbours[1] != -1) {
        std::vector<double> send = data.get_col(domain.imax);
        std::vector<double> recv;
        MPI_Sendrecv(&send, domain.size_y + 2, MPI_DOUBLE, domain.neighbours[1], 2, &recv, domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[1], 1, MPI_COMM_WORLD, &status);
        data.set_col(recv, domain.imax + 1);
    }

    // Communicate to Top
    if (domain.neighbours[2] != -1) {
        std::vector<double> send = data.get_row(domain.jmax);
        std::vector<double> recv;
        MPI_Sendrecv(&send, domain.size_x + 2, MPI_DOUBLE, domain.neighbours[2], 3, &recv, domain.size_x + 2,
                     MPI_DOUBLE, domain.neighbours[2], 4, MPI_COMM_WORLD, &status);
        data.set_row(recv, domain.jmax + 1);
    }

    // Communicate to Bottom
    if (domain.neighbours[3] != -1) {
        std::vector<double> send = data.get_row(1);
        std::vector<double> recv;
        MPI_Sendrecv(&send, domain.size_x + 2, MPI_DOUBLE, domain.neighbours[3], 4, &recv, domain.size_x + 2,
                     MPI_DOUBLE, domain.neighbours[3], 3, MPI_COMM_WORLD, &status);
        data.set_row(recv, 0);
    }
}

double reduce_min(double dt) {
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double final_dt;

    MPI_Reduce(&dt, &final_dt, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(&final_dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return final_dt;
}