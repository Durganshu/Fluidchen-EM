#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "Datastructures.hpp"
#include "Domain.hpp"
#include "Enums.hpp"

/**
 * @brief Data structure holds cells and related sub-containers
 *
 */
class Grid {
  public:
    Grid() = default;

    /**
     * @brief Constructor for the Grid
     *
     * @param[in] geometry file name
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] cell size in x direction
     * @param[in] cell size in y direction
     *
     */
    Grid(std::string geom_name, Domain &domain, int iproc, int jproc);

    /// index based cell access
    Cell cell(int i, int j) const;

    /// access number of cells in x direction
    int imax() const;
    /// access number of cells in y direction
    int jmax() const;
    /// access number of cells in x direction including ghost cells
    int imaxb() const;
    /// access number of cells in x direction including ghost cells
    int jmaxb() const;

    /// access number of cells in x direction excluding ghost cells
    const Domain &domain() const;

    /// access cell size in x-direction
    double dx() const;
    /// access cell size in y-direction
    double dy() const;

    /**
     * @brief Access fluid cells
     *
     * @param[out] vector of fluid cells
     */
    const std::vector<Cell *> &fluid_cells() const;

    /**
     * @brief Access inflow cells
     *
     * @param[out] vector of fluid cells
     */
    const std::vector<Cell *> &inflow_cells() const;

    /**
     * @brief Access outflow cells
     *
     * @param[out] vector of fluid cells
     */
    const std::vector<Cell *> &outflow_cells() const;


    /**
     * @brief Access moving wall cells
     *
     * @param[out] vector of moving wall cells
     */
    const std::vector<Cell *> &moving_wall_cells() const;

    /**
     * @brief Access fixed wall cells
     *
     * @param[out] vector of fixed wall cells
     */
    const std::vector<Cell *> &fixed_wall_cells() const;

    /**
     * @brief Access cold fixed wall cells
     *
     * @param[out] vector of cold fixed wall cells
     */
    const std::vector<Cell *> &cold_fixed_wall_cells() const;

    /**
     * @brief Access hot fixed wall cells
     *
     * @param[out] vector of hot fixed wall cells
     */
    const std::vector<Cell *> &hot_fixed_wall_cells() const;

    /**
     * @brief Access adiabatic fixed wall cells
     *
     * @param[out] vector of adiabatic fixed wall cells
     */
    const std::vector<Cell *> &adiabatic_fixed_wall_cells() const;

  private:
    /**@brief Default lid driven cavity case generator
     *
     * This function creates default lid driven cavity
     * case without need for a pgm file
     */
    void build_lid_driven_cavity();

    /// Build cell data structures with given geometrical data
    void assign_cell_types(std::vector<std::vector<int>> &geometry_data);
    /// Extract geometry from pgm file and create geometrical data
    void parse_geometry_file(std::string filedoc, std::vector<std::vector<int>> &geometry_data);

    Matrix<Cell> _cells;
    std::vector<Cell *> _fluid_cells;
    std::vector<Cell *> _inflow_cells;
    std::vector<Cell *> _outflow_cells;
    std::vector<Cell *> _fixed_wall_cells;
    std::vector<Cell *> _cold_fixed_wall_cells;
    std::vector<Cell *> _hot_fixed_wall_cells;
    std::vector<Cell *> _adiabatic_fixed_wall_cells;
    std::vector<Cell *> _moving_wall_cells;

    Domain _domain;

    double _dx;
    double _dy;

    int _rank;
    int _size;
    int _iproc;
    int _jproc;
};
