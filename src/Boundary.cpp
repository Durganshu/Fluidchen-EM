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
            break;
        }

        else if (j == 0) {
            field.u(i, j) = -field.u(i, j + 1);
            field.v(i, j) = 0;
            break;
        }

        else if (i == field.p_matrix().imax()) {
            field.u(i, j) = 0;
            field.v(i + 1, j) = field.v(i, j);
            break;
        }

        else if (j == field.p_matrix().jmax()) {
            field.u(i, j + 1) = -field.u(i, j);
            field.v(i, j) = 0;
            break;
        }
        else std::cout<<"Error in FixedWallBoundary::apply() "<<"/n";
        //if (i==field.)
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {}
