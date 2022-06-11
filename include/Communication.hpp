#include <mpi.h>
#include "Datastructures.hpp"
#include "Domain.hpp"

/**
 * @brief Communicates ghost layer values with neighbours
 *
 * @param[in] Field to be applied
 */
static void communicate(Matrix<double> &data, Domain &domain);

/**
 * @brief Returns the minimum time-step among all processes
 *
 * @param[in] time-steps
 */
double reduce_min(double dt);