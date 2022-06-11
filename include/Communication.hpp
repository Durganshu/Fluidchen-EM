#include <mpi.h>
#include "Datastructures.hpp"

/**
 * @brief Communicates ghost layer values with neighbours
 *
 * @param[in] Field to be applied
 */
static void communicate(Matrix &data, Domain &domain);

/**
 * @brief Returns the minimum time-step among all processes
 *
 * @param[in] time-steps
 */
double reduce_min(double dt);