#include "Discretization.hpp"
#include <cmath>
#include <iostream>
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

// Convection in x direction
double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {

    double du2dx;
    double duvdy;

    du2dx = (U(i, j) + U(i + 1, j)) * (U(i, j) + U(i + 1, j)) - 
            (U(i - 1, j) + U(i, j)) * (U(i - 1, j) + U(i, j));
    du2dx += _gamma * (std::abs(U(i, j) + U(i + 1, j)) * (U(i, j) - U(i + 1, j)) -
                       std::abs(U(i - 1, j) + U(i, j)) * ((U(i - 1, j) - U(i, j))));
    du2dx *= (0.25 / _dx);

    duvdy =
        (V(i, j) + V(i + 1, j)) * (U(i, j) + U(i, j + 1)) - (V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) + U(i, j));
    duvdy += _gamma * (std::abs(V(i, j) + V(i + 1, j)) * (U(i, j) - U(i, j + 1)) -
                       std::abs(V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) - U(i, j)));
    duvdy *= (0.25 / _dy);

    return (du2dx + duvdy);
}

// Convection in y direction
double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {

    double duvdx;
    double dv2dy;

    duvdx = (U(i, j) + U(i, j + 1)) * (V(i, j) + V(i + 1, j)) - (U(i - 1, j) + 
            U(i - 1, j + 1)) * (V(i - 1, j) + V(i, j));
    duvdx += _gamma * (std::abs(U(i, j) + U(i, j + 1)) * (V(i, j) - V(i + 1, j))
         - std::abs(U(i - 1, j) + U(i - 1, j + 1)) * (V(i - 1, j) - V(i, j)));
    duvdx *= (0.25 / _dx);

    dv2dy = (V(i, j) + V(i, j + 1)) * (V(i, j) + V(i, j + 1)) - 
            (V(i, j - 1) + V(i, j)) * (V(i, j - 1) + V(i, j));
    dv2dy += _gamma * (std::abs(V(i, j) + V(i, j + 1)) * (V(i, j) - V(i, j + 1))
        - std::abs(V(i, j - 1) + V(i, j)) * (V(i, j - 1) - V(i, j)));
    dv2dy *= (0.25 / _dy);

    return (duvdx + dv2dy);
}

double Discretization::convection_t(const Matrix<double> &U, const Matrix<double> &V, const Matrix<double> &T, int i,
                                    int j) {

    double duTdx = U(i, j) * (T(i, j) + T(i + 1, j)) - U(i - 1, j) * (T(i - 1, j) + T(i, j));
    duTdx += _gamma * (std::abs(U(i, j)) * (T(i, j) - T(i + 1, j)) - 
        std::abs(U(i - 1, j)) * (T(i - 1, j) - T(i, j)));
    duTdx *= (0.5 / _dx);

    double dvTdy = V(i, j) * (T(i, j) + T(i, j + 1)) - V(i, j - 1) * (T(i, j - 1) + T(i, j));
    dvTdy += _gamma * (std::abs(V(i, j)) * (T(i, j) - T(i, j + 1)) - 
            std::abs(V(i, j - 1)) * (T(i, j - 1) - T(i, j)));
    dvTdy *= (0.5 / _dy);

    return (duTdx + dvTdy);
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}
