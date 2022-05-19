#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>

Fields::Fields(Grid &grid, double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {

    _U = Matrix<double>(imax + 2, jmax + 2);
    _V = Matrix<double>(imax + 2, jmax + 2);
    _P = Matrix<double>(imax + 2, jmax + 2);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();

        _U(i, j) = UI;
        _V(i, j) = VI;
        _P(i, j) = PI;
    }
}

void Fields::calculate_fluxes(Grid &grid) {
    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();
        if (i < grid.imax())
            _F(i, j) = _gx + _U(i, j) +
                       _dt * ((_nu * Discretization::laplacian(_U, i, j)) - Discretization::convection_u(_U, _V, i, j));
        if (j < grid.jmax())
            _G(i, j) = _gy + _V(i, j) +
                       _dt * ((_nu * Discretization::laplacian(_V, i, j)) - Discretization::convection_v(_U, _V, i, j));
    }
}

void Fields::calculate_rs(Grid &grid) {
    auto idt = 1. / _dt; // Calculate 1/dt
    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();
        rs(i, j) = idt * (((_F(i, j) - _F(i - 1, j)) / grid.dx()) + ((_G(i, j) - _G(i, j - 1)) / grid.dy()));
    }
}

void Fields::calculate_velocities(Grid &grid) {

    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();

        if (i < grid.imax()) _U(i, j) = _F(i, j) - (_dt / grid.dx()) * (_P(i + 1, j) - _P(i, j));

        if (j < grid.jmax()) _V(i, j) = _G(i, j) - (_dt / grid.dy()) * (_P(i, j + 1) - _P(i, j));
    }
}

double Fields::calculate_dt(Grid &grid) {

    auto max_u = 0.0;
    auto max_v = 0.0;

    for (auto &elem : grid.fluid_cells()) {

        int i = elem->i();
        int j = elem->j();

        if (std::fabs(_U(i, j)) > max_u) max_u = std::fabs(_U(i, j));

        if (std::fabs(_V(i, j)) > max_v) max_v = std::fabs(_V(i, j));
    }

    auto factor1 = 1 / (1 / (grid.dx() * grid.dx()) + 1 / (grid.dy() * grid.dy()));
    factor1 = factor1 / (2 * _nu);
    auto factor2 = grid.dx() / max_u;
    auto factor3 = grid.dy() / max_v;
    _dt = _tau * std::min(factor1, std::min(factor2, factor3));
    return _dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
