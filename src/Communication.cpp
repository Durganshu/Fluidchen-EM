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
    //static int count = 0;  RANK =0

    // std::ofstream rank_0, rank_1;
    
    // rank_0.open("rank_0.csv");
    // rank_1.open("rank_1.csv");
    // Communicate to Left
    //MPI_Barrier(MPI_COMM_WORLD);
    if (domain.neighbours[0] != -1) {
        // std::cout << "In left " << rank << " \n";
        std::vector<double> send = data.get_col(1);
        std::vector<double> recv(domain.size_y + 2);

        MPI_Status status;
        
        MPI_Sendrecv(send.data(), domain.size_y + 2, MPI_DOUBLE, domain.neighbours[0], 1, &recv[0], domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[0], 2, MPI_COMM_WORLD, &status);
        

        data.set_col(recv, 0);

        
        // std::cout << "Printing receive (right) for rank = "<< rank << ", count = "<< count<< " \n";
        // for (auto & elem : recv){
        //     if(rank == 0) rank_0 << "(" <<elem << ", " << rank << ")\n";
        //     else rank_1 << "(" <<elem << ", " << rank << ")\n";
        //     //std::cout << elem << "\n";
        // }
    }

    // Communicate to Right
    //MPI_Barrier(MPI_COMM_WORLD);
    if (domain.neighbours[1] != -1) {

        // std::cout << "In right " << rank << " " << domain.imax << " \n";
        std::vector<double> send = data.get_col(domain.size_x);
        std::vector<double> recv(domain.size_y + 2);

        
        // std::cout << "Printing send (left) for rank = "<< rank << ", count = "<< count<< " \n";
        // for (auto & elem : send){
        //     if(rank == 0) rank_0 << "(" <<elem << ", " << rank << ")\n";
        //     else rank_1 << "(" <<elem << ", " << rank << ")\n";
        //     //std::cout << elem << "\n";
        // }
        MPI_Status status;
        MPI_Sendrecv(send.data(), domain.size_y + 2, MPI_DOUBLE, domain.neighbours[1], 2, &recv[0], domain.size_y + 2,
                     MPI_DOUBLE, domain.neighbours[1], 1, MPI_COMM_WORLD, &status);
        data.set_col(recv, domain.size_x + 1);

        
    }
    // Communicate to Top
    if (domain.neighbours[2] != -1) {
        // std::cout << "In top " << rank << " \n";
        std::vector<double> send = data.get_row(domain.size_y);
        std::vector<double> recv(domain.size_x + 2);
        MPI_Status status;
        MPI_Sendrecv(send.data(), domain.size_x + 2, MPI_DOUBLE, domain.neighbours[2], 3, &recv[0], domain.size_x + 2,
                     MPI_DOUBLE, domain.neighbours[2], 4, MPI_COMM_WORLD, &status);
        data.set_row(recv, domain.size_y + 1);
    }

    // Communicate to Bottom
    if (domain.neighbours[3] != -1) {
        // std::cout << "In bottom " << rank << " \n";
        std::vector<double> send = data.get_row(1);
        std::vector<double> recv(domain.size_x + 2);
        MPI_Status status;
        MPI_Sendrecv(send.data(), domain.size_x + 2, MPI_DOUBLE, domain.neighbours[3], 4, &recv[0], domain.size_x + 2,
                     MPI_DOUBLE, domain.neighbours[3], 3, MPI_COMM_WORLD, &status);
        data.set_row(recv, 0);
    }
    //count++;
}


double Communication::reduce_min(double x)
{


    // double final_dt;

    // MPI_Reduce(&dt,&final_dt,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    // MPI_Bcast(&final_dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    // return final_dt;

    double final_x;
    MPI_Allreduce(&x, &final_x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return final_x;


}

double Communication::reduce_sum(double res) {
    double reduced_res;
    ///REPLACE WITH ALL_REDUCE LIKEABOVE
    MPI_Allreduce(&res, &reduced_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Bcast(&reduced_res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return reduced_res;
}
