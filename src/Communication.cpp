#include "Communication.hpp"
#include <mpi.h>

void Communication::finalize()
{
    MPI_Finalize();
}

void Communication::init_parallel(int *argn, char **args, int &rank, int &size)
{
    MPI_Init(argn,&args);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

void Communication::communicate(Matrix<double> &data, const Domain &domain, int rank) {

    //Communicate to Left
    if (domain.neighbours[0] != -1) {
        std::vector<double> send = data.get_col(1);
        std::vector<double> recv(domain.size_y + 2);
        MPI_Status status;
        MPI_Sendrecv(send.data(), domain.size_y + 2, MPI_DOUBLE, domain.neighbours[0], 1, &recv[0], domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[0], 2, MPI_COMM_WORLD, &status);
        data.set_col(recv, 0);
    }
    // Communicate to Right
    if (domain.neighbours[1] != -1) {
        std::vector<double> send = data.get_col(domain.size_x);
        std::vector<double> recv(domain.size_y + 2);
        MPI_Status status;
        MPI_Sendrecv(send.data(), domain.size_y + 2, MPI_DOUBLE, domain.neighbours[1], 2, &recv[0], domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[1], 1, MPI_COMM_WORLD, &status);
        data.set_col(recv, domain.size_x + 1);   
    }
    // Communicate to Top
    if (domain.neighbours[2] != -1) {
        std::vector<double> send = data.get_row(domain.size_y);
        std::vector<double> recv(domain.size_x + 2);
        MPI_Status status;
        MPI_Sendrecv(send.data(), domain.size_x + 2, MPI_DOUBLE, domain.neighbours[2], 3, &recv[0], domain.size_x + 2,
                     MPI_DOUBLE, domain.neighbours[2], 4, MPI_COMM_WORLD, &status);
        data.set_row(recv, domain.size_y + 1);
    }
    // Communicate to Bottom
    if (domain.neighbours[3] != -1) {;
        std::vector<double> send = data.get_row(1);
        std::vector<double> recv(domain.size_x + 2);
        MPI_Status status;
        MPI_Sendrecv(send.data(), domain.size_x + 2, MPI_DOUBLE, domain.neighbours[3], 4, &recv[0], domain.size_x + 2,
                     MPI_DOUBLE, domain.neighbours[3], 3, MPI_COMM_WORLD, &status);
        data.set_row(recv, 0);
    }
}


double Communication::reduce_min(double x)
{
    double final_x;
    MPI_Allreduce(&x, &final_x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return final_x;
}

double Communication::reduce_sum(double res) {
    double reduced_res;
    MPI_Allreduce(&res, &reduced_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return reduced_res;
}
