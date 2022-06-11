#include "Datastructures.hpp"
#include "Domain.hpp"
#include <mpi.h>


/**
 * @brief Communicates ghost layer values with neighbours
 *
 * @param[in] Field to be applied
 */
void communicate(Matrix<double> &data, const Domain &domain);

/**
 * @brief Returns the minimum time-step among all processes
 *
 * @param[in] time-steps
 */
double reduce_min(double dt);