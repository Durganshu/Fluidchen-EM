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
        //std::cout<<"i,j="<<i<<","<<j<<"\n";
        
        // TOP implies that the top border of the cell exists i.e.
        // these cells should be in the "bottommost row"
        if (elem->is_border(border_position::TOP)) {
            field.u(i, j) = -field.u(i, j + 1);
            field.v(i, j) = 0;
            field.p(i, j) = field.p(i , j + 1);
            field.g(i, j) = field.v(i, j);
            //std::cout<<"i,j="<<i<<","<<j<<"\n";
        }

        // RIGHT implies that the right border of the cell exists i.e.
        // these cells should be in the "leftmost column"
        else if (elem->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = -field.v(i + 1, j);
            field.p(i, j) = field.p(i + 1, j);
            field.f(i, j) = field.u(i, j);
            //std::cout<<"i,j="<<i<<","<<j<<"\n";
        }

        // LEFT implies that the left border of the cell exists i.e.
        // these cells should be in the "rightmost column"
        else if (elem->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i - 1, j) = -field.v(i, j);
            field.p(i - 1, j) = field.p(i, j);
            field.f(i - 1, j) = field.u(i - 1, j);
            //std::cout<<"i,j="<<i<<","<<j<<"\n";
        }

        //There won't be any top row in the fixed wall boundary
        // else if (elem->is_border(border_position::BOTTOM)) {
        //     field.u(i, j + 1) = -field.u(i, j);
        //     field.v(i, j) = 0;
        //     field.p(i, j + 1) = field.p(i, j);
        //     std::cout<<"i,j="<<i<<","<<j<<"\n";
        // } 
        else{
                      //<< "\n"; // Throw error if the index is not on boundary
            // std::cout << "Corner cells:";
            // std::cout<<"i,j="<<i<<","<<j<<"\n";
            }//<< "\n"; // Throw error if the index is not on boundary

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
        //std::cout<<"i,j="<<i<<","<<j<<"\n";
        //field.u(i, j) = 2*1 - field.u(i, j - 1);
        //std::cout<<_wall_velocity.begin()->second<<"\n";
        field.u(i, j) = 2 * (_wall_velocity.begin()->second) - field.u(i, j - 1);
        field.v(i, j - 1) = 0;
        field.p(i, j) = field.p(i, j - 1);
        field.g(i, j - 1) = field.v(i, j - 1);
        //std::cout<<"i,j="<<i<<","<<j<<" : u = "<< field.u(i, j)<<"\n";
    }
}
