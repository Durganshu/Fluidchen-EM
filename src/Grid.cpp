#include "Grid.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

Grid::Grid(std::string geom_name, Domain &domain, int iproc, int jproc) {


    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    _domain = domain;
    _iproc= iproc;
    _jproc= jproc;

    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

    if (geom_name.compare("NONE")) {
        std::vector<std::vector<int>> geometry_data(_domain.size_x + 2,
                                                    std::vector<int>(_domain.size_y + 2, 0));
        parse_geometry_file(geom_name, geometry_data);
        assign_cell_types(geometry_data);
    } else {
        build_lid_driven_cavity();
    }
}

void Grid::build_lid_driven_cavity() {
    std::vector<std::vector<int>> geometry_data(_domain.size_x + 2,
                                                std::vector<int>(_domain.size_y + 2, 0));

    for (int i = 0; i < _domain.size_x + 2; ++i) {
        for (int j = 0; j < _domain.size_y + 2; ++j) {
            // Bottom, left and right walls: no-slip
            if ((i == 0 && _domain.imin ==0)||(i==_domain.size_x +1 && _domain.imax==_domain.domain_size_x +2 )|| (j == 0  && _domain.jmin==0) ) {
                geometry_data.at(i).at(j) = LidDrivenCavity::fixed_wall_id;
            }
            // Top wall: moving wall
            else if (j == _domain.size_y + 1 && _domain.jmax==_domain.domain_size_y+2) {
                geometry_data.at(i).at(j) = LidDrivenCavity::moving_wall_id;
            }
        }
    }
    // if (_rank == 0) {   //uncomment to visualize the part of the domain for rank 0
    //     std::cout<< std::endl;
    //     std::cout<<_domain.size_x<<"  "<<_domain.size_y<<std::endl;
    //     for (int j = _domain.size_y + 1; j >= 0; --j) {
    //         for (int i = 0; i < _domain.size_x + 2; ++i) {
    //             std::cout << geometry_data.at(i).at(j) << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }
    
    assign_cell_types(geometry_data);
}

void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data) {

    int i = 0;
    int j = 0;

    for (int j_geom = 0; j_geom < _domain.size_y+2; ++j_geom) { //modified limits to account _cells for each process
        { i = 0; }
        for (int i_geom = 0; i_geom < _domain.size_x+2; ++i_geom) {//modified limits to account _cells for each process
            if (geometry_data.at(i_geom).at(j_geom) == 0) {
                _cells(i, j) = Cell(i, j, cell_type::FLUID);
                _fluid_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 1) {
                _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data.at(i_geom).at(j_geom));
                _inflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 2) {
                _cells(i, j) = Cell(i, j, cell_type::OUTFLOW, geometry_data.at(i_geom).at(j_geom));
                _outflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 3) {
                _cells(i, j) = Cell(i, j, cell_type::COLD_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _cold_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 4) {
                _cells(i, j) = Cell(i, j, cell_type::HOT_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _hot_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 5) {
                _cells(i, j) = Cell(i, j, cell_type::ADIABATIC_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _adiabatic_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 6) {
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == LidDrivenCavity::moving_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data.at(i_geom).at(j_geom));
                _moving_wall_cells.push_back(&_cells(i, j));
            }

            ++i;
        }
        ++j;
    }

    // Corner cell neighbour assigment
    // Bottom-Left Corner
    i = 0;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
    }
    // Top-Left Corner
    i = 0;
    j = _domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
    }

    // Top-Right Corner
    i = _domain.size_x + 1;
    j = Grid::_domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
    }
    // Bottom-Right Corner
    i = Grid::_domain.size_x + 1;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
    }
    // Bottom cells
    j = 0;
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).type() != cell_type::FLUID) {
            if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::RIGHT);
            }
            if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::LEFT);
            }
            if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::TOP);
            }
        }
    }
    // Top Cells
    j = Grid::_domain.size_y + 1;

    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        if (_cells(i, j).type() != cell_type::FLUID) {
            if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::RIGHT);
            }
            if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::LEFT);
            }
            if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::BOTTOM);
            }
        }
    }

    // Left Cells
    i = 0;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).type() != cell_type::FLUID) {
            if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::RIGHT);
            }
            if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::BOTTOM);
            }
            if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::TOP);
            }
        }
    }
    // Right Cells
    i = Grid::_domain.size_x + 1;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);

        if (_cells(i, j).type() != cell_type::FLUID) {
            if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::LEFT);
            }
            if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::BOTTOM);
            }
            if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::TOP);
            }
        }
    }

    // Inner cells
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        for (int j = 1; j < _domain.size_y + 1; ++j) {
            _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
            _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
            _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
            _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);

            if (_cells(i, j).type() != cell_type::FLUID) {
                if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::LEFT);
                }
                if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::RIGHT);
                }
                if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::BOTTOM);
                }
                if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::TOP);
                }
            }
        }
    }

}

void Grid::parse_geometry_file(std::string filedoc, std::vector<std::vector<int>> &geometry_data) {
     if (_rank == 0) {

    int numcols, numrows, depth;
    //Vector to store whole geometry which is read by only rank 0
    std::vector<std::vector<int>> entire_geometry_data(_domain.domain_size_x + 2, std::vector<int>(_domain.domain_size_y + 2, 0));
    std::ifstream infile(filedoc);
    std::stringstream ss;
    std::string inputLine = "";

    // First line : version
    getline(infile, inputLine);
    if (inputLine.compare("P2") != 0) {
        std::cerr << "First line of the PGM file should be P2" << std::endl;
    }

    // Second line : comment
    getline(infile, inputLine);

    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numrows >> numcols;
    // Fourth line : depth
    ss >> depth;

    // Following lines : data
    for (int col = numcols - 1; col > -1; --col) {
        for (int row = 0; row < numrows; ++row) {
            ss >> entire_geometry_data[row][col];
        }
    }

    infile.close();
    

    int I, J;
    int imin, jmin, imax, jmax;

    // Sending Data to other processors
    for (int i = 1; i < _size; ++i) {
        I = i % _iproc + 1;
        J = i / _iproc + 1;
        imin = (I - 1) * ((numrows - 2) / _iproc);
        imax = I * ((numrows - 2) / _iproc) + 2;
        jmin = (J - 1) * ((numcols - 2) / _jproc);
        jmax = J * ((numcols - 2) / _jproc) + 2;

        // Adding the extra cells when number of cells is not divisible by iproc and jproc
        if (I == _iproc) imax = numrows;

        if (J == _jproc) jmax = numcols;


        std::vector<int> rank_geometry_data;
        for (int row = imin; row < imax; ++row) {
            for (int col = jmin; col < jmax; ++col) {
                rank_geometry_data.push_back(entire_geometry_data[row][col]);
            }
        }
        // std::cout<<"\n product "<<((imax-imin)*(jmax-jmin))<<std::endl;
        // //Send to each rank
        // std::cout<<"Sent Size "<<rank_geometry_data.size()<<std::endl;
        MPI_Send(rank_geometry_data.data(), rank_geometry_data.size(), MPI_INT, i, 999999, MPI_COMM_WORLD);
    }

        // Assigning Geometry data for rank 0
    for (int col = 0; col < (_domain.size_y + 2); ++col) {
        for (int row = 0; row < (_domain.size_x + 2); ++row) {
            geometry_data[row][col] = entire_geometry_data[row][col];
        }
    }


} //End of if where rank 0 is working
//************************************************************************************************************
/// Replace _domain.imax - _domain.imin by size_x+2 for uniformity everywhere
else {
        // Receive data from rank 0
        std::vector<int> rank_geometry_data((_domain.size_x + 2) * (_domain.size_y + 2),0);
        std::cout<<"Receive Size"<<rank_geometry_data.size()<<std::endl;
        MPI_Status status;

        MPI_Recv(rank_geometry_data.data(), rank_geometry_data.size(), MPI_INT, 0, 999999, MPI_COMM_WORLD, &status);

        std::cout << " Size of geometry_data: " << geometry_data.size() << ", rank = "<< _rank << "\n";
        std::cout << " Size of rank_geometry_data: " << rank_geometry_data.size() << ", rank = "<< _rank << " \n";
        for (int col = 0; col < _domain.size_y + 2; ++col) {
            for (int row = 0; row < _domain.size_x + 2; ++row) {
                geometry_data.at(row).at(col) = rank_geometry_data[row * (_domain.size_y + 2) + col];
            }
        }
    }



}// End of Parse Geometry

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

double Grid::dx() const { return _domain.dx; }

double Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

const std::vector<Cell *> &Grid::fluid_cells() const { return _fluid_cells; }

const std::vector<Cell *> &Grid::inflow_cells() const { return _inflow_cells; }

const std::vector<Cell *> &Grid::outflow_cells() const { return _outflow_cells; }

const std::vector<Cell *> &Grid::fixed_wall_cells() const { return _fixed_wall_cells; }

const std::vector<Cell *> &Grid::moving_wall_cells() const { return _moving_wall_cells; }

const std::vector<Cell *> &Grid::cold_fixed_wall_cells() const { return _cold_fixed_wall_cells; }

const std::vector<Cell *> &Grid::hot_fixed_wall_cells() const { return _hot_fixed_wall_cells; }

const std::vector<Cell *> &Grid::adiabatic_fixed_wall_cells() const { return _adiabatic_fixed_wall_cells; }
