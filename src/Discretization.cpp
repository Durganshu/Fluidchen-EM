#include "Discretization.hpp"
#include <iostream>
#include <cmath>
#include <math.h>
#include <type_traits>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}
//Convection in x direction
double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // auto val = (1 / _dx) * ((((U(i, j) + U(i + 1, j)) * (U(i, j) + U(i + 1, j))) * 0.25) -
    //                     (((U(i - 1, j) + U(i, j)) * (U(i - 1, j) + U(i, j))) * 0.25)) +
    //        _gamma * (1 / _dx) *
    //            ((((std::fabs(U(i, j) + U(i + 1, j)) * ((U(i, j) - U(i + 1, j)))) * 0.25)) -
    //             ((std::fabs(U(i - 1, j) + U(i, j)) * (U(i - 1, j) - U(i, j))) * 0.25)) +
    //        (1 / _dy) * (((((V(i, j) + V(i + 1, j)) * (U(i, j) + U(i, j + 1))) * 0.25) -
    //                      ((V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) + U(i, j))) * 0.25)) +
    //        _gamma * (1 / _dy) *
    //            ((((std::fabs(V(i, j) + V(i + 1, j)) * ((U(i, j) - U(i, j + 1)))) * 0.25)) -
    //             ((std::fabs(V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) - U(i, j))) * 0.25));


    //std::cout<<val<<"\n";
    return (1 / _dx) * ((((U(i, j) + U(i + 1, j)) * (U(i, j) + U(i + 1, j))) * 0.25) -
                        (((U(i - 1, j) + U(i, j)) * (U(i - 1, j) + U(i, j))) * 0.25)) +
           _gamma * (1 / _dx) *
               ((((std::fabs(U(i, j) + U(i + 1, j)) * ((U(i, j) - U(i + 1, j)))) * 0.25)) -
                ((std::fabs(U(i - 1, j) + U(i, j)) * (U(i - 1, j) - U(i, j))) * 0.25)) +
           (1 / _dy) * (((((V(i, j) + V(i + 1, j)) * (U(i, j) + U(i, j + 1))) * 0.25) -
                         ((V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) + U(i, j))) * 0.25)) +
           _gamma * (1 / _dy) *
               ((((std::fabs(V(i, j) + V(i + 1, j)) * ((U(i, j) - U(i, j + 1)))) * 0.25)) -
                ((std::fabs(V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) - U(i, j))) * 0.25));
}
//Convection in y direction
double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {

    return (1 / _dx) * ((((U(i, j) + U(i , j+1)) * (V(i, j) + V(i + 1, j))) * 0.25) -
                        (((U(i - 1, j) + U(i-1, j+1)) * (V(i - 1, j) + V(i, j))) * 0.25)) +
           _gamma * (1 / _dx) *
               ((((std::fabs(U(i, j) + U(i , j+1)) * ((V(i, j) - V(i + 1, j)))) * 0.25)) -
                ((std::fabs(U(i - 1, j) + U(i-1, j+1)) * (V(i - 1, j) - V(i, j))) * 0.25)) +
           (1 / _dy) * (((((V(i, j) + V(i , j+1)) * (V(i, j) + V(i, j + 1))) * 0.25) -
                         ((V(i, j - 1) + V(i , j)) * (V(i, j - 1) + V(i, j))) * 0.25)) +
           _gamma * (1 / _dy) *
               ((((std::fabs(V(i, j) + V(i , j+1)) * ((V(i, j) - V(i, j + 1)))) * 0.25)) -
                ((std::fabs(V(i, j - 1) + V(i , j )) * (V(i, j - 1) - V(i, j))) * 0.25));
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {}
