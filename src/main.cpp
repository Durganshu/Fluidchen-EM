#include <iostream>
#include <string>
#include <omp.h>

#include "Case.hpp"

void printIntro();

int main(int argn, char **args) {

    if (argn > 1) {
        int rank,size;
        MPI_Init(&argn,&args);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        #define MPI_RANK rank
        #define MPI_SIZE size
        std::string file_name{args[1]};
        Case problem(file_name, argn, args, rank);
        problem.printIntro();
        problem.simulate();
        MPI_Finalize();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}