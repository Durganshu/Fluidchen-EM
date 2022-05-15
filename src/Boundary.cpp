#include "Boundary.hpp"
#include "Enums.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {

    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();

        // TOP implies that the top border of the cell exists i.e.
        // these cells should be in the "bottommost row"
        if (elem->is_border(border_position::TOP)) {
            field.u(i, j) = -field.u(i, j + 1);
            field.v(i, j) = 0;
            field.p(i, j) = field.p(i, j + 1);
            field.g(i, j) = field.v(i, j);
        }

        // RIGHT implies that the right border of the cell exists i.e.
        // these cells should be in the "leftmost column"
        else if (elem->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = -field.v(i + 1, j);
            field.p(i, j) = field.p(i + 1, j);
            field.f(i, j) = field.u(i, j);
        }

        // LEFT implies that the left border of the cell exists i.e.
        // these cells should be in the "rightmost column"
        else if (elem->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = -field.v(i - 1, j);
            field.p(i, j) = field.p(i - 1, j);
            field.f(i - 1, j) = field.u(i - 1, j);
        }
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
        field.u(i, j) = 2 * (_wall_velocity.begin()->second) - field.u(i, j - 1);
        field.v(i, j - 1) = 0;
        field.p(i, j) = field.p(i, j - 1);
        field.g(i, j - 1) = field.v(i, j - 1);
    }
}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double inflow_x_velocity, double inflow_y_velocity)
    : _cells(cells),_x_velocity(inflow_x_velocity),_y_velocity(inflow_y_velocity) {}
void InflowBoundary::apply(Fields &field) {}
OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells, double outflow_pressure)
    : _cells(cells),_outflow_pressure(outflow_pressure) {}
void OutflowBoundary::apply(Fields &field) {}