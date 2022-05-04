#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {

    // for (auto &elem : _cells) {
    //     int i = elem.i();
    //     int j = elem.j();

    //     if (i == 0 || i == _U.imax())) {
    //         if (j == 0) {
    //             _U[i][j] = -_U[i][j + 1];
    //             _V[i][j] = -_V[i+1][j]
    //             return;
    //         }
    //         if (j == _U.jmax()) {
    //             _U[i][j] = -_U[i][j - 1];
                
    //         }conda 
    //     }
    // }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {}
