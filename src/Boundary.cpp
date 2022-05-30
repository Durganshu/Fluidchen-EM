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

        if (check_neighbours(elem) > 2) {
            std::cout << "Boundary cell at i = " << i << ", j = " << j
                      << " has more than two fluid cells as neighbours. Please fix the geometry file. Exiting!"
                      << "\n";
            exit(0);
        }

        if (elem->is_border(border_position::TOP)) {

            // NE corner
            if (elem->is_border(border_position::RIGHT)) {
                field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.u(elem->neighbour(border_position::LEFT)->i(), j) =
                    -field.u(elem->neighbour(border_position::LEFT)->i(), elem->neighbour(border_position::TOP)->j());
                field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = -field.v(
                    elem->neighbour(border_position::RIGHT)->i(), elem->neighbour(border_position::BOTTOM)->j());
            }

            // NW corner
            else if (elem->is_border(border_position::LEFT)) {

                field.u(elem->neighbour(border_position::LEFT)->i(), j) = 0;
                field.v(i, j) = 0;
                field.u(i, j) = -field.u(i, elem->neighbour(border_position::TOP)->j());
                field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = -field.v(
                    elem->neighbour(border_position::LEFT)->i(), elem->neighbour(border_position::BOTTOM)->j());

            }

            // Special case : A cell having both TOP and BOTTOM boundaries
            else if (elem->is_border(border_position::BOTTOM)) {
                field.u(i, j) = -field.u(i, elem->neighbour(border_position::BOTTOM)->j());
                field.v(i, j) = 0;
                field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = 0;
                field.u(elem->neighbour(border_position::LEFT)->i(), j) =
                    -(field.u(elem->neighbour(border_position::LEFT)->i(),
                              elem->neighbour(border_position::BOTTOM)->j()) +
                      field.u(elem->neighbour(border_position::LEFT)->i(),
                              elem->neighbour(border_position::TOP)->j())) /
                    2;
            }

            // Cells only having TOP boundary
            else {

                field.u(i, j) = -field.u(i, elem->neighbour(border_position::TOP)->j());
                field.v(i, j) = 0;
            }

        }

        else if (elem->is_border(border_position::BOTTOM)) {

            // SE corner
            if (elem->is_border(border_position::RIGHT)) {
                field.u(i, j) = 0;
                field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = 0;
                field.u(elem->neighbour(border_position::LEFT)->i(), j) = -field.u(
                    elem->neighbour(border_position::LEFT)->i(), elem->neighbour(border_position::BOTTOM)->j());
                field.v(i, j) = -field.v(elem->neighbour(border_position::RIGHT)->i(), j);
                2;
            }

            // SW corner
            else if (elem->is_border(border_position::LEFT)) {
                field.u(elem->neighbour(border_position::LEFT)->i(), j) = 0;
                field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = 0;
                field.u(i, j) = -field.u(i, elem->neighbour(border_position::BOTTOM)->j());
                field.v(i, j) = -field.v(elem->neighbour(border_position::LEFT)->i(), j);
            }

            // Cells only having BOTTOM boundary
            else {
                field.u(i, j) = -field.u(i, elem->neighbour(border_position::BOTTOM)->j());
                field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = 0;
            }
        }

        else if (elem->is_border(border_position::RIGHT)) {

            // Special case : A cell having both RIGHT and LEFT boundaries
            if (elem->is_border(border_position::LEFT)) {
                field.u(i, j) = 0;
                field.u(elem->neighbour(border_position::LEFT)->i(), j) = 0;
                field.v(i, j) = -(field.v(elem->neighbour(border_position::RIGHT)->i(), j) +
                                  field.v(elem->neighbour(border_position::LEFT)->i(), j)) /
                                2;
            }

            // Cells only having RIGHT boundary
            else {
                field.u(i, j) = 0;
                field.v(i, j) = -field.v(elem->neighbour(border_position::RIGHT)->i(), j);
            }
        }

        // Cells only having LEFT boundary
        else if (elem->is_border(border_position::LEFT)) {
            field.u(elem->neighbour(border_position::LEFT)->i(), j) = 0;
            field.v(i, j) = -field.v(elem->neighbour(border_position::LEFT)->i(), j);
        }
    }
}

void FixedWallBoundary::apply_pressure(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();

        if (elem->is_border(border_position::TOP)) {

            if (elem->is_border(border_position::RIGHT)) {
                field.p(i, j) = (field.p(i, elem->neighbour(border_position::TOP)->j()) +
                                 field.p(elem->neighbour(border_position::RIGHT)->i(), j)) /
                                2;
            } else if (elem->is_border(border_position::LEFT)) {
                field.p(i, j) = (field.p(i, elem->neighbour(border_position::TOP)->j()) +
                                 field.p(elem->neighbour(border_position::LEFT)->i(), j)) /
                                2;
            } else if (elem->is_border(border_position::BOTTOM)) { // Need to verify
                field.p(i, j) = (field.p(i, elem->neighbour(border_position::TOP)->j()) +
                                 field.p(i, elem->neighbour(border_position::BOTTOM)->j())) /
                                2;
            } else {
                field.p(i, j) = field.p(i, elem->neighbour(border_position::TOP)->j());
            }

        }

        else if (elem->is_border(border_position::BOTTOM)) {

            if (elem->is_border(border_position::RIGHT)) {
                field.p(i, j) = (field.p(elem->neighbour(border_position::RIGHT)->i(), j) +
                                 field.p(i, elem->neighbour(border_position::BOTTOM)->j())) /
                                2;
            } else if (elem->is_border(border_position::LEFT)) {
                field.p(i, j) = (field.p(i, elem->neighbour(border_position::BOTTOM)->j()) +
                                 field.p(elem->neighbour(border_position::LEFT)->i(), j)) /
                                2;
            } else {
                field.p(i, j) = field.p(i, elem->neighbour(border_position::BOTTOM)->j());
            }
        }

        else if (elem->is_border(border_position::RIGHT)) {

            if (elem->is_border(border_position::LEFT)) {
                field.p(i, j) = (field.p(elem->neighbour(border_position::RIGHT)->i(), j) +
                                 field.p(elem->neighbour(border_position::LEFT)->i(), j)) /
                                2;
            } else {
                field.p(i, j) = field.p(elem->neighbour(border_position::RIGHT)->i(), j);
            }
        }

        else if (elem->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(elem->neighbour(border_position::LEFT)->i(), j);
        }
    }
}

void FixedWallBoundary::apply_temperature(Fields &field) const {

    // Storing the value of wall temperature
    const double wall_id = _wall_temperature.begin()->first;
    const double wall_temperature = _wall_temperature.begin()->second;

    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();

        // Applying the temperature BCs according to the wall id
        if (elem->is_border(border_position::TOP)) {
            // For NE corner cell
            if (elem->is_border(border_position::RIGHT)) {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) = 2 * wall_temperature - (field.t(i, elem->neighbour(border_position::TOP)->j()) +
                                                            field.t(elem->neighbour(border_position::RIGHT)->i(), j)) *
                                                               0.5;
                } else {

                    field.t(i, j) = 0.5 * (field.t(i, elem->neighbour(border_position::TOP)->j()) +
                                           field.t(elem->neighbour(border_position::RIGHT)->i(), j));
                }
            }

            // For NW corner cell (Fluid is present in the north and west of this corner cell)
            else if (elem->is_border(border_position::LEFT)) {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) = 2 * wall_temperature - (field.t(i, elem->neighbour(border_position::TOP)->j()) +
                                                            field.t(elem->neighbour(border_position::LEFT)->i(), j)) *
                                                               0.5;
                } else {

                    field.t(i, j) = 0.5 * (field.t(i, elem->neighbour(border_position::TOP)->j()) +
                                           field.t(elem->neighbour(border_position::LEFT)->i(), j));
                }
            }

            //  (Fluid is present in both north and south of this cell)
            else if (elem->is_border(border_position::BOTTOM)) {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) = 2 * wall_temperature - (field.t(i, elem->neighbour(border_position::TOP)->j()) +
                                                            field.t(i, elem->neighbour(border_position::BOTTOM)->j())) *
                                                               0.5;
                } else {

                    field.t(i, j) = 0.5 * (field.t(i, elem->neighbour(border_position::TOP)->j()) +
                                           field.t(i, elem->neighbour(border_position::BOTTOM)->j()));
                }
            }

            // For bottommost boundary  (Fluid is present ONLY  in the north  of these cells)
            else {
                if (wall_id == 3 || wall_id == 4) {

                    field.t(i, j) = 2 * wall_temperature - field.t(i, elem->neighbour(border_position::TOP)->j());
                } else {
                    field.t(i, j) = field.t(i, elem->neighbour(border_position::TOP)->j());
                }
            }
        }

        if (elem->is_border(border_position::BOTTOM)) {

            // SE cell
            if (elem->is_border(border_position::RIGHT)) {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) = 2 * wall_temperature - (field.t(i, elem->neighbour(border_position::BOTTOM)->j()) +
                                                            field.t(elem->neighbour(border_position::RIGHT)->i(), j)) *
                                                               0.5;
                } else {

                    field.t(i, j) = 0.5 * (field.t(i, elem->neighbour(border_position::BOTTOM)->j()) +
                                           field.t(elem->neighbour(border_position::RIGHT)->i(), j));
                }
            }

            // For SW corner cell (Fluid is present in the north and west of this corner cell)
            else if (elem->is_border(border_position::LEFT)) {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) = 2 * wall_temperature - (field.t(i, elem->neighbour(border_position::BOTTOM)->j()) +
                                                            field.t(elem->neighbour(border_position::LEFT)->i(), j)) *
                                                               0.5;
                } else {

                    field.t(i, j) = 0.5 * (field.t(i, elem->neighbour(border_position::BOTTOM)->j()) +
                                           field.t(elem->neighbour(border_position::LEFT)->i(), j));
                }
            }

            // For topmost boundary
            else {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) = 2 * wall_temperature - field.t(i, elem->neighbour(border_position::BOTTOM)->j());
                } else {
                    field.t(i, j) = field.t(i, elem->neighbour(border_position::BOTTOM)->j());
                }
            }
        }

        if (elem->is_border(border_position::LEFT)) {

            // Fluid cells exist on the left and right borders
            if (elem->is_border(border_position::RIGHT)) {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) =
                        2 * wall_temperature - 0.5 * (field.t(elem->neighbour(border_position::LEFT)->i(), j) +
                                                      field.t(elem->neighbour(border_position::RIGHT)->i(), j));
                } else {
                    field.t(i, j) = 0.5 * (field.t(elem->neighbour(border_position::LEFT)->i(), j) +
                                           field.t(elem->neighbour(border_position::RIGHT)->i(), j));
                }

            }
            // For rightmost boundary (Fluid cells exist on the left)
            else {
                if (wall_id == 3 || wall_id == 4) {
                    field.t(i, j) = 2 * wall_temperature - field.t(elem->neighbour(border_position::LEFT)->i(), j);
                } else {
                    field.t(i, j) = field.t(elem->neighbour(border_position::LEFT)->i(), j);
                }
            }
        }
        // For leftmost boundary
        if (elem->is_border(border_position::RIGHT)) {
            if (wall_id == 3 || wall_id == 4) {
                field.t(i, j) = 2 * wall_temperature - field.t(elem->neighbour(border_position::RIGHT)->i(), j);
            }

            else {
                field.t(i, j) = field.t(elem->neighbour(border_position::RIGHT)->i(), j);
            }
        }
    }
}

int FixedWallBoundary::check_neighbours(Cell *cell) {
    int number_of_fluid_neighbours = 0;
    if (cell->is_border(border_position::TOP) && cell->neighbour(border_position::TOP)->type() == cell_type::FLUID)
        number_of_fluid_neighbours++;
    if (cell->is_border(border_position::BOTTOM) &&
        cell->neighbour(border_position::BOTTOM)->type() == cell_type::FLUID)
        number_of_fluid_neighbours++;
    if (cell->is_border(border_position::LEFT) && cell->neighbour(border_position::LEFT)->type() == cell_type::FLUID)
        number_of_fluid_neighbours++;
    if (cell->is_border(border_position::RIGHT) && cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID)
        number_of_fluid_neighbours++;

    return number_of_fluid_neighbours;
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
        field.u(i, j) =
            2 * (_wall_velocity.begin()->second) - field.u(i, elem->neighbour(border_position::BOTTOM)->j());
        field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = 0;
    }
}

void MovingWallBoundary::apply_pressure(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();
        field.p(i, j) = field.p(i, elem->neighbour(border_position::BOTTOM)->j());
    }
}

// Temperature BC for moving_wall is not in the scope of this worksheet so the function is kept as a dummy.
void MovingWallBoundary::apply_temperature(Fields &field) const {}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double inflow_x_velocity, double inflow_y_velocity)
    : _cells(cells), _x_velocity(inflow_x_velocity), _y_velocity(inflow_y_velocity) {}

void InflowBoundary::apply(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();

        field.u(i, j) = _x_velocity;
        field.v(i, j) = -field.v(elem->neighbour(border_position::RIGHT)->i(), j);
    }
}

void InflowBoundary::apply_pressure(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();
        field.p(i, j) = field.p(elem->neighbour(border_position::RIGHT)->i(), j);
    }
}

// Temperature BC for inflow is not in the scope of this worksheet so the function is kept as a dummy.
void InflowBoundary::apply_temperature(Fields &field) const {}

OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells, double outlet_pressure) : _cells(cells), _pressure(outlet_pressure) {}

void OutflowBoundary::apply(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();

        field.u(i, j) = field.u(elem->neighbour(border_position::LEFT)->i(), j);
        field.v(i, j) = field.v(elem->neighbour(border_position::LEFT)->i(), j);
    }
}

void OutflowBoundary::apply_pressure(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();
        field.p(i, j) = _pressure;
    }
}

// Temperature BC for outflow is not in the scope of this worksheet so the function is kept as a dummy.
void OutflowBoundary::apply_temperature(Fields &field) const {}