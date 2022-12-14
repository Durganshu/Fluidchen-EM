#include "Grid.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

Grid::Grid(std::string geom_name, Domain &domain, int iproc, int jproc, int rank, int size) {

    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    _domain = domain;
    _iproc = iproc;
    _jproc = jproc;
    _rank = rank;
    _size = size;

    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

    if (geom_name.compare("NONE")) {
        std::vector<std::vector<int>> geometry_data(_domain.size_x + 2, std::vector<int>(_domain.size_y + 2, 0));
        parse_geometry_file(geom_name, geometry_data);
        assign_cell_types(geometry_data);
    } else {
        build_lid_driven_cavity();
    }
}

void Grid::build_lid_driven_cavity() {
    std::vector<std::vector<int>> geometry_data(_domain.size_x + 2, std::vector<int>(_domain.size_y + 2, 0));

    for (int i = 0; i < _domain.size_x + 2; ++i) {
        for (int j = 0; j < _domain.size_y + 2; ++j) {
            // Bottom, left and right walls: no-slip
            if ((i == 0 && _domain.imin == 0) ||
                (i == _domain.size_x + 1 && _domain.imax == _domain.domain_size_x + 2) ||
                (j == 0 && _domain.jmin == 0)) {
                geometry_data.at(i).at(j) = LidDrivenCavity::fixed_wall_id;
            }
            // Top wall: moving wall
            else if (j == _domain.size_y + 1 && _domain.jmax == _domain.domain_size_y + 2) {
                geometry_data.at(i).at(j) = LidDrivenCavity::moving_wall_id;
            }
        }
    }
    assign_cell_types(geometry_data);
}

void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data) {

    int i = 0;
    int j = 0;

    for (int j_geom = 0; j_geom < _domain.size_y + 2; ++j_geom) {
        { i = 0; }
        for (int i_geom = 0; i_geom < _domain.size_x + 2; ++i_geom) {

            // Set to true if the cell is a buffer cell
            bool isBuffer = false;

            if (i_geom == 0 && _domain.neighbours[0] != -1) {
                isBuffer = true;
            } else if (i_geom == (_domain.size_x + 1) && _domain.neighbours[1] != -1) {
                isBuffer = true;
            } else if (j_geom == (_domain.size_y + 1) && _domain.neighbours[2] != -1) {
                isBuffer = true;
            } else if (j_geom == 0 && _domain.neighbours[3] != -1) {
                isBuffer = true;
            }

            if (geometry_data.at(i_geom).at(j_geom) == 0) {
                _cells(i, j) = Cell(i, j, cell_type::FLUID);
                if (!isBuffer) _fluid_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 1) {
                _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _inflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 2) {
                _cells(i, j) = Cell(i, j, cell_type::OUTFLOW, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _outflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 3) {
                _cells(i, j) = Cell(i, j, cell_type::COLD_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _cold_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 4) {
                _cells(i, j) = Cell(i, j, cell_type::HOT_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _hot_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 5) {
                _cells(i, j) = Cell(i, j, cell_type::ADIABATIC_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _adiabatic_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 6) {
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == LidDrivenCavity::moving_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _moving_wall_cells.push_back(&_cells(i, j));

            //********Assign higher potential cell type
            } else if (geometry_data.at(i_geom).at(j_geom) == 7) {
                _cells(i, j) = Cell(i, j, cell_type::HIGHER_POTENTIAL_WALL, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _higher_potential_cells.push_back(&_cells(i, j));
                if (!isBuffer) _fixed_wall_cells.push_back(&_cells(i, j));
            }  //*********Assign lower potential cell type
            else if (geometry_data.at(i_geom).at(j_geom) == 8) {
                _cells(i, j) = Cell(i, j, cell_type::LOWER_POTENTIAL_WALL, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _lower_potential_cells.push_back(&_cells(i, j));
                if (!isBuffer) _fixed_wall_cells.push_back(&_cells(i, j));
            }
            else if(geometry_data.at(i_geom).at(j_geom) == 9){
                _cells(i, j) = Cell(i, j, cell_type::FLUID, geometry_data.at(i_geom).at(j_geom));
                if (!isBuffer) _coupled_cells.push_back(&_cells(i, j));
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
        // Vector to store whole geometry which is read by only rank 0
        std::vector<std::vector<int>> entire_geometry_data(_domain.domain_size_x + 2,
                                                           std::vector<int>(_domain.domain_size_y + 2, 0));
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

            // Adding/Removing the extra cells when number of cells is not divisible by iproc and jproc
            if (I == _iproc) imax = numrows;
            if (J == _jproc) jmax = numcols;

            // Assign geometry data for each rank
            std::vector<int> rank_geometry_data;
            for (int row = imin; row < imax; ++row) {
                for (int col = jmin; col < jmax; ++col) {
                    rank_geometry_data.push_back(entire_geometry_data[row][col]);
                }
            }

            // Send geometry data to each rank
            MPI_Send(rank_geometry_data.data(), rank_geometry_data.size(), MPI_INT, i, 999999, MPI_COMM_WORLD);
        }

        // Assigning Geometry data for rank 0
        for (int col = 0; col < (_domain.size_y + 2); ++col) {
            for (int row = 0; row < (_domain.size_x + 2); ++row) {
                geometry_data[row][col] = entire_geometry_data[row][col];
            }
        }
    }
    else {
        // Receive data from rank 0
        std::vector<int> rank_geometry_data((_domain.size_x + 2) * (_domain.size_y + 2), 0);
        
        MPI_Status status;

        MPI_Recv(rank_geometry_data.data(), rank_geometry_data.size(), MPI_INT, 0, 999999, MPI_COMM_WORLD, &status);

        for (int col = 0; col < _domain.size_y + 2; ++col) {
            for (int row = 0; row < _domain.size_x + 2; ++row) {
                geometry_data.at(row).at(col) = rank_geometry_data[row * (_domain.size_y + 2) + col];
            }
        }
    }

} // End of Parse Geometry

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

const std::vector<Cell *> &Grid::higher_potential_cells() const { return _higher_potential_cells; }

const std::vector<Cell *> &Grid::lower_potential_cells() const { return _lower_potential_cells; }

const std::vector<Cell *> &Grid::coupled_cells() const { return _coupled_cells; }
