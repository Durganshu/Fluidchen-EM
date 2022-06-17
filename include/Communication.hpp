#include <mpi.h>
#include "Fields.hpp"

class Communication {
  public:
    /**
     * @brief Constructor for the Communication.
     */
    Communication() = default;
     /**
     * @brief Initialize MPI Parallelization
     *
     *
     * @param[in] Input file name
     */
    static void init_parallel(int* argn,char** args,int &rank, int &size);
    /**
     * @brief Finalize MPI 
     */
    static void finalize();
     /**
     * @brief Function to communicate field data
     * Communicates field data from one domain to another
     * @param[in] Field data to be communicated
     * @param[in] Domain
     */
    static void communicate(Matrix<double> &data, const Domain &domain, int rank);
    /**
     * @brief Function to find the minimum value among all
     * Finds the minimum among all domain and broadcasts it to all
     * @param[in] local dt of each rank
     */
    static double reduce_min(double dt);
    /**
     * @brief Function to sum a variable across all domains
     * 
     * @param[in] local residual of each rank
     */
    static double reduce_sum(double res);

  
};