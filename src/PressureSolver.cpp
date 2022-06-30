#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.imax() + 1 && j != grid.jmax() + 1) { // exclude the buffer cells

            field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                            coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
        }
    }

    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.imax() + 1 && j != grid.jmax() + 1) { // exclude the buffer cells
            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }
    }
    // Every domain just returns the rloc
    return rloc;
}

double SOR::solve_potential(Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundary) {

    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

    for (auto &elem : boundary) {
        elem->apply_potential(phi_matrix());
    }

    // Run an SOR Iteration
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.imax() + 1 && j != grid.jmax() + 1) { // exclude the buffer cells
            phi(i, j) = (1.0 - _omega) * phi(i, j) + coeff * (Discretization::sor_helper(phi_matrix(), i, j));
        }
    }

    // Calcualte Residual
    double res = 0.0;
    double rloc = 0.0;
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.imax() + 1 && j != grid.jmax() + 1) { // exclude the buffer cells
            double val = Discretization::laplacian(phi_matrix(), i, j);
            rloc += (val * val);
        }
    }
    return rloc = std::sqrt(rloc / grid.fluid_cells().size()); // Final residual
}