#include "Fields.hpp"

#include <algorithm>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {}

void Fields::calculate_rs(Grid &grid) {
    auto idt = 1. / _dt; // Calculate 1/dt
    for (auto i = 1; i <= grid.imax(); i++) {
        for (auto j = 1; j <= grid.jmax(); j++) {
            rs(i,j) = (idt) *((f(i,j)-f(i-1,j)/grid.dx()+(g(i,j)-g(i,j-1)/grid.dy())));
        }
    }
}

void Fields::calculate_velocities(Grid &grid) {
    for (auto i = 1; i <= grid.imax(); i++) {
        for (auto j = 1; j <= grid.jmax(); j++) {
            u(i, j) = f(i, j) - (_dt / grid.dx()) * (p(i + 1, j) - p(i, j));
            v(i, j) = g(i, j) - (_dt / grid.dy()) * (p(i, j + 1) - p(i, j));
        }
    }
}

double Fields::calculate_dt(Grid &grid) { 
    
    auto max_u = 0.0;
    auto max_v = 0.0;

    for (auto &elem : grid.fluid_cells()) {
    
        int i = elem->i();
        int j = elem->j();

        if(u(i, j) > max_u) max_u = u(i, j);

        if(v(i, j) > max_v) max_v = v(i, j);
    }
    
    auto factor1 = 1/(1/(grid.dx()*grid.dx()) + 1/(grid.dy()*grid.dy()));
    factor1 = factor1/(2 * _nu);
    auto factor2 = grid.dx()/max_u;
    auto factor3 = grid.dy()/max_v;
    _dt = _tau*std::min(factor1, std::min(factor2, factor3));
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

