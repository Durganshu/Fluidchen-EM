#include<mpi.h>

#include"Fields.hpp"

/**
 * @brief 
 * 
 * @param[in] Field to be applied
 */
 void communicate(Matrix<double> &data, const Domain &domain);

 double reduce_min(double dt);

 double reduce_sum(double res);