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
 * @brief Returns the minimum parameter across all processes
 *
 * @param[in] Parameter to be reduced to minimum
 */
double reduce_min(double x);

/**
 * @brief Returns the maximum parameter across all processes
 *
 * @param[in] Parameter to be reduced to maximum
 */
double reduce_max(double x);

/**
 * @brief Returns the total of given parameter across all processes
 *
 * @param[in] Parameter to be reduced to total
 */
double reduce_total(double x);

/**
 * @brief Returns the total of given parameter across all processes
 *
 * @param[in] Parameter to be reduced to total
 */
int reduce_total(int x);