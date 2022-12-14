#include <iostream>
#include <string>
#include <chrono>

#include "Case.hpp"
#include "Communication.hpp"

//void printIntro();

int main(int argn, char **args) {

    auto start = std::chrono::steady_clock::now();

    int rank, size;
    Communication::init_parallel(&argn, args, rank, size);

    if (argn > 1) {
        if (rank == 0) Case::printIntro();
        std::string file_name{args[1]};
        Case problem(file_name, argn, args, rank, size);
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
    if(rank==0){
    std::cout << "\nSimulation Complete!\n";
    auto end = std::chrono::steady_clock::now();
    std::cout << "Software Runtime:" << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s\n\n";
    }
    

    return 0;
}