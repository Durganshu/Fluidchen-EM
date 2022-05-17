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
        
        if (check_neighbours(elem) > 2){
            std::cout<<"Boundary cell at i = "<< i<<", j = "<<j<<
                " has more than two fluid cells as neighbours. Please fix the geometry file. Exiting!"<< "\n";
            exit(0);
        }
        
        // TOP implies that the top border of the cell exists i.e.
        // these cells should be in the "bottommost row"
        if (elem->is_border(border_position::TOP)) {
            //std::cout << "i = " << i<<", "<<"j = "<<j<<"\n";
             //std::cout << "i = " << i<<", "<<"j = "<<elem->neighbour(border_position::TOP)->j()<<"\n";
            field.u(i, j) = -field.u(i, elem->neighbour(border_position::TOP)->j());
            field.v(i, j) = 0;
            field.p(i, j) = field.p(i, elem->neighbour(border_position::TOP)->j());
            field.g(i, j) = field.v(i, j);
        }

        // BOTTOM implies that the bottom border of the cell exists i.e.
        // these cells should be in the "topmost row"
       if (elem->is_border(border_position::BOTTOM)) {
            //std::cout << "i = " << i<<", "<<"j = "<<elem->neighbour(border_position::BOTTOM)->j()<<"\n";
            field.u(i, j) = -field.u(i, elem->neighbour(border_position::BOTTOM)->j());
            field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = 0;
            field.p(i, j) = field.p(i, elem->neighbour(border_position::BOTTOM)->j());
            field.g(i, elem->neighbour(border_position::BOTTOM)->j()) = field.v(i, elem->neighbour(border_position::BOTTOM)->j());
        } 

        // RIGHT implies that the right border of the cell exists i.e.
        // these cells should be in the "leftmost column"
        else if (elem->is_border(border_position::RIGHT)) {
            //std::cout << "i = " << i<<", "<<"j = "<<j<<"\n";
            field.u(i, j) = 0;
            field.v(i, j) = -field.v(elem->neighbour(border_position::RIGHT)->i(), j);
            field.p(i, j) = field.p(elem->neighbour(border_position::RIGHT)->i(), j);
            field.f(i, j) = field.u(i, j);
        }

        // LEFT implies that the left border of the cell exists i.e.
        // these cells should be in the "rightmost column"
        else if (elem->is_border(border_position::LEFT)) {
            //std::cout << "i = " << elem->neighbour(border_position::LEFT)->i()<<", "<<"j = "<<j<<"\n";
            field.u(elem->neighbour(border_position::LEFT)->i(), j) = 0;
            field.v(i, j) = -field.v(elem->neighbour(border_position::LEFT)->i(), j);
            field.p(i, j) = field.p(elem->neighbour(border_position::LEFT)->i(), j);
            field.f(elem->neighbour(border_position::LEFT)->i(), j) = field.u(elem->neighbour(border_position::LEFT)->i(), j);
        }
    }
}
int FixedWallBoundary::check_neighbours(Cell * cell){
    int number_of_fluid_neighbours = 0;
    if (cell->is_border(border_position::TOP) && 
        cell->neighbour(border_position::TOP)->type() == cell_type::FLUID)
        number_of_fluid_neighbours++;
    if (cell->is_border(border_position::BOTTOM) && 
        cell->neighbour(border_position::BOTTOM)->type() == cell_type::FLUID)
        number_of_fluid_neighbours++;
    if (cell->is_border(border_position::LEFT) && 
        cell->neighbour(border_position::LEFT)->type() == cell_type::FLUID)
        number_of_fluid_neighbours++;
    if (cell->is_border(border_position::RIGHT) && 
        cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID)
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
        //std::cout << "i = " << i<<", "<<"j = "<<j<<"\n";
        //std::cout << "i = " << i<<", "<<"j = "<<elem->neighbour(border_position::BOTTOM)->j()<<"\n";
        field.u(i, j) = 2 * (_wall_velocity.begin()->second) - field.u(i, elem->neighbour(border_position::BOTTOM)->j());
        field.v(i, elem->neighbour(border_position::BOTTOM)->j()) = 0;
        field.p(i, j) = field.p(i, elem->neighbour(border_position::BOTTOM)->j());
        //field.g(i, j - 1) = field.v(i, elem->neighbour(border_position::BOTTOM)->j());
    }
}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double inflow_x_velocity, double inflow_y_velocity, double inflow_pressure)
    : _cells(cells),_x_velocity(inflow_x_velocity),_y_velocity(inflow_y_velocity), _pressure(inflow_pressure) {}

void InflowBoundary::apply(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();
        
        field.u(i, j) = _x_velocity;
        field.v(i, j) = 2*_y_velocity-field.v(elem->neighbour(border_position::RIGHT)->i(),j);
        field.p(i, j) = _pressure;
        field.f(i, j) = field.u(i, j);
        //field.g(i, j) = field.v(i, j);
        
    }
}

OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells)
    : _cells(cells) {}

void OutflowBoundary::apply(Fields &field) {
    for (auto &elem : _cells) {
        int i = elem->i();
        int j = elem->j();
        
        field.u(elem->neighbour(border_position::LEFT)->i(), j)=field.u(i, j) ; 
        field.v(i, j) = field.v(elem->neighbour(border_position::LEFT)->i(), j);
        field.p(i, j) = 0.0; 
        field.f(elem->neighbour(border_position::LEFT)->i(), j) = field.u(elem->neighbour(border_position::LEFT)->i(), j);
        //field.g(i, j) = field.v(i, j);
       
    
    }
}