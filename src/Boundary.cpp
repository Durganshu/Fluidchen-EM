#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {

    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();

        if (i == 0) {
            field.u(i, j) = 0;
            field.v(i, j) = -field.v(i + 1, j);
            field.p(i, j) = field.p(i + 1, j);
            break;
        }

        else if (j == 0) {
            field.u(i, j) = -field.u(i, j + 1);
            field.v(i, j) = 0;
            field.p(i, j) = field.p(i, j + 1);
            break;
        }

        else if (i == field.p_matrix().imax()) {
            field.u(i, j) = 0;
            field.v(i + 1, j) = field.v(i, j);
            field.p(i + 1, j) = field.p(i, j);
            break;
        }

        else if (j == field.p_matrix().jmax()) {
            field.u(i, j + 1) = -field.u(i, j);
            field.v(i, j) = 0;
            field.p(i, j + 1) = field.p(i, j);
            break;
        } else 
            std::cout << "Error in FixedWallBoundary::apply() "
                      << "/n"; //Throw error if the index is on boundary
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();
        field.u(i, j) = 1;
        field.v(i, j) = 0;

        //Check if the y-index is jmax
        if (j != field.p_matrix().jmax())
            std::cout << "Error in MovingWallBoundary::apply() "
                      << "/n";
    }
}
