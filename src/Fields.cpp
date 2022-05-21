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

Fields::Fields(Grid &grid, double nu, double alpha, double beta, double dt, double tau, int imax, int jmax, double UI,
               double VI, double PI, double TI)
    : _nu(nu), _alpha(alpha), _beta(beta), _dt(dt), _tau(tau) {

    _U = Matrix<double>(imax + 2, jmax + 2);
    _V = Matrix<double>(imax + 2, jmax + 2);
    _P = Matrix<double>(imax + 2, jmax + 2);
    _T = Matrix<double>(imax + 2, jmax + 2);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();

        _U(i, j) = UI;
        _V(i, j) = VI;
        _P(i, j) = PI;
        _T(i, j) = TI;
    }
}

void Fields::calculate_temperatures(Grid &grid){};

void Fields::calculate_fluxes(Grid &grid) {
    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();
        //if (i < grid.imax())
            _F(i, j) =  _U(i, j) +
                       _dt * ((_nu * Discretization::laplacian(_U, i, j)) - Discretization::convection_u(_U, _V, i, j) + _gx);
        //if (j < grid.jmax())
            _G(i, j) = _V(i, j) +
                       _dt * ((_nu * Discretization::laplacian(_V, i, j)) - Discretization::convection_v(_U, _V, i, j) + _gy);
    }

    ////Applying Flux BC to fixed walls
    for (auto &elem : grid.fixed_wall_cells()) {
        int i = elem->i();
        int j = elem->j();

     if (elem->is_border(border_position::TOP)) {
            // std::cout << "i = " << i<<", "<<"j = "<<j<<"\n";
            // std::cout << "i = " << i<<", "<<"j = "<<elem->neighbour(border_position::TOP)->j()<<"\n";
            if (elem->is_border(border_position::RIGHT)) {
                _F(i, j) = 0.0;
                _G(i, j) = 0.0;
            }

            else if (elem->is_border(border_position::LEFT)) {
                _F(elem->neighbour(border_position::LEFT)->i(), j) = 0.0;
                //_F(i - 1, j) = 0.0;
                _G(i, j) = 0.0;

            }

            else if (elem->is_border(border_position::BOTTOM)) { // Need to verify
          
                _G(i, j) = 0.0;
                //_G(i, j - 1) = 0.0;
                _G(i, elem->neighbour(border_position::BOTTOM)->j()) = 0.0;

            }

            else {
              
                _G(i, j) = _V(i, j);
            }

        }

        // BOTTOM implies that the bottom border of the cell exists i.e.
        // these cells should be in the "topmost row"
        else if (elem->is_border(border_position::BOTTOM)) {
            // std::cout << "i = " << i<<", "<<"j = "<<elem->neighbour(border_position::BOTTOM)->j()<<"\n";
            if (elem->is_border(border_position::RIGHT)) {

                _F(i, j) = 0.0;
                //_G(i, j - 1) = 0.0;
                _G(i, elem->neighbour(border_position::BOTTOM)->j()) = 0.0;
            }


            else if (elem->is_border(border_position::LEFT)) {
 
                //_F(i - 1, j) = 0.0;
                //_G(i, j - 1) = 0.0;
                _F(elem->neighbour(border_position::LEFT)->i(), j) = 0.0;
                _G(i, elem->neighbour(border_position::BOTTOM)->j()) = 0.0;
            }

            else {
              
                _G(i, elem->neighbour(border_position::BOTTOM)->j()) =
                _V(i, elem->neighbour(border_position::BOTTOM)->j());

                //_G(i,j-1) =_V(i, j-1);
            }
        }

        // RIGHT implies that the right border of the cell exists i.e.
        // these cells should be in the "leftmost column"
        else if (elem->is_border(border_position::RIGHT)) {
            // std::cout << "i = " << i<<", "<<"j = "<<j<<"\n";
            if (elem->is_border(border_position::LEFT)) { // Need to verify
               
                _F(i, j) = 0.0;
                //_F(i - 1, j) = 0.0;
                _F(elem->neighbour(border_position::LEFT)->i(), j) = 0.0;

            }

            else {
                
                _F(i, j) = _U(i, j);
            }
        }

        // LEFT implies that the left border of the cell exists i.e.
        // these cells should be in the "rightmost column"
        else if (elem->is_border(border_position::LEFT)) {
            // std::cout << "i = " << elem->neighbour(border_position::LEFT)->i()<<", "<<"j = "<<j<<"\n";

            _F(elem->neighbour(border_position::LEFT)->i(), j) =
                 _U(elem->neighbour(border_position::LEFT)->i(), j);

            //_F(i-1, j) =_U(i-1, j);
        }
    
    
    }

    // Flux setup for Moving wall cells
    for (auto &elem : grid.moving_wall_cells()) {
        int i = elem->i();
        int j = elem->j();

        _G(i, elem->neighbour(border_position::BOTTOM)->j()) = _V(i, elem->neighbour(border_position::BOTTOM)->j());
        //_G(i, j - 1) = _V(i, elem->neighbour(border_position::BOTTOM)->j());
    }

    // Flux setup for inflow cells
    for (auto &elem : grid.inflow_cells()) {
        int i = elem->i();
        int j = elem->j();

        _F(i, j) = _U(i, j);
        //_G(i, j) = _V(i, j);  Not required
    }

    // Flux setup for outflow cells
    for (auto &elem : grid.outflow_cells()) {
        int i = elem->i();
        int j = elem->j();

        _F(elem->neighbour(border_position::LEFT)->i(), j) =
            _U(elem->neighbour(border_position::LEFT)->i(), j);
        //_F(i-1, j) =_U(i-1, j);
    }

}

void Fields::calculate_rs(Grid &grid) {
    auto idt = 1. / _dt; // Calculate 1/dt
    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();
        rs(i, j) = idt * (((_F(i, j) - _F(elem->neighbour(border_position::LEFT)->i(), j)) / grid.dx()) + 
                                        ((_G(i, j) - _G(i, elem->neighbour(border_position::BOTTOM)->j())) / grid.dy()));
        //rs(i, j) = idt * (((_F(i, j) - _F(i - 1, j)) / grid.dx()) + ((_G(i, j) - _G(i, j - 1)) / grid.dy()));
    }
}

void Fields::calculate_velocities(Grid &grid) {

    for (const auto &elem : grid.fluid_cells()) {
        int i = elem->i();
        int j = elem->j();

        _U(i, j) = _F(i, j) - (_dt / grid.dx()) * (_P(elem->neighbour(border_position::RIGHT)->i(), j) - _P(i, j)); 
        //_U(i, j) = _F(i, j) - (_dt / grid.dx()) * (_P(i + 1, j) - _P(i, j));

        _V(i, j) = _G(i, j) - (_dt / grid.dy()) * (_P(i, elem->neighbour(border_position::TOP)->j()) - _P(i, j));
        //_V(i, j) = _G(i, j) - (_dt / grid.dy()) * (_P(i, j + 1) - _P(i, j));
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

double Fields::calculate_dt_e(Grid &grid) {

    auto max_u = 0.0;
    auto max_v = 0.0;

    for (auto &elem : grid.fluid_cells()) {

        int i = elem->i();
        int j = elem->j();

        if (std::fabs(_U(i, j)) > max_u) max_u = std::fabs(_U(i, j));

        if (std::fabs(_V(i, j)) > max_v) max_v = std::fabs(_V(i, j));
    }

    auto factor = 1 / (1 / (grid.dx() * grid.dx()) + 1 / (grid.dy() * grid.dy()));
    auto factor1 = factor / (2 * _nu);
    auto factor2 = grid.dx() / max_u;
    auto factor3 = grid.dy() / max_v;
    auto factor4 = factor/(2*_alpha);
    _dt = _tau * std::min(std::min(factor1, std::min(factor2, factor3)),factor4);
    return _dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::t(int i, int j) { return _T(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
