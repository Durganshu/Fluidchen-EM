#include "Communication.hpp"
#include <vector>

void communicate(Matrix<double> &data, const Domain &domain) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    // Communicate to Left
    if (domain.neighbours[0] != -1) {
        std::cout << "In left " << rank << " \n";
        std::vector<double> send = data.get_col(1);
        std::vector<double> recv(domain.size_y + 2);

        MPI_Sendrecv(send.data(), domain.size_y + 2, MPI_DOUBLE, domain.neighbours[0], 1, &recv[0], domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[0], 2, MPI_COMM_WORLD, &status);

        data.set_col(recv, 0);
    }

    // Communicate to Right
    if (domain.neighbours[1] != -1) {
        std::cout << "In right " << rank << " " << domain.imax << " \n";
        std::vector<double> send = data.get_col(domain.size_x);
        std::vector<double> recv(domain.size_y + 2);
        MPI_Sendrecv(send.data(), domain.size_y + 2, MPI_DOUBLE, domain.neighbours[1], 2, &recv[0], domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[1], 1, MPI_COMM_WORLD, &status);
        data.set_col(recv, domain.size_x + 1);
    }

    // Communicate to Top
    if (domain.neighbours[2] != -1) {
        std::cout << "In top " << rank << " \n";
        std::vector<double> send = data.get_row(domain.size_y);
        std::vector<double> recv(domain.size_x + 2);
        MPI_Sendrecv(send.data(), domain.size_x + 2, MPI_DOUBLE, domain.neighbours[2], 3, &recv[0], domain.size_x + 2,
                     MPI_DOUBLE, domain.neighbours[2], 4, MPI_COMM_WORLD, &status);
        data.set_row(recv, domain.size_y + 1);
    }

    // Communicate to Bottom
    if (domain.neighbours[3] != -1) {
        std::cout << "In bottom " << rank << " \n";
        std::vector<double> send = data.get_row(1);
        std::vector<double> recv(domain.size_x + 2);
        MPI_Sendrecv(send.data(), domain.size_x + 2, MPI_DOUBLE, domain.neighbours[3], 4, &recv[0], domain.size_x + 2,
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