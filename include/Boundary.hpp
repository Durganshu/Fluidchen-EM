#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"

/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Main method to patch the boundary conditons to given field and
     * grid
     *
     * @param[in] Field to be applied
     */
    virtual void apply(Fields &field) = 0;
    virtual void apply_pressure(Fields &field) = 0;
    virtual void apply_temperature(Fields &field) = 0;
    virtual ~Boundary() = default;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature);
    int check_neighbours(Cell *cell);
    virtual ~FixedWallBoundary() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressure(Fields &field);
    void apply_temperature(Fields &field);
    bool _temperature_boundary_type;   ////****************NOT REQUIRED ANYMORE
  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_temperature;
};

/**
 * @brief Inflow boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities and pressure
 */
class InflowBoundary : public Boundary {
  public:
    InflowBoundary(std::vector<Cell *> cells, double inflow_x_velocity, double inflow_y_velocity,
                   double inflow_pressure);
    virtual ~InflowBoundary() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressure(Fields &field);
    bool _temperature_boundary_type;  ////****************NOT REQUIRED ANYMORE
    virtual void apply_temperature(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _x_velocity;
    double _y_velocity;
    double _pressure;
};

/**
 * @brief Outflow boundary condition for the outer boundaries of the domain.
 * Dirichlet for pressure, Neumann for velocities
 */
class OutflowBoundary : public Boundary {
  public:
    OutflowBoundary(std::vector<Cell *> cells);
    virtual ~OutflowBoundary() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressure(Fields &field);
    bool _temperature_boundary_type ;  ////****************NOT REQUIRED ANYMORE
    virtual void apply_temperature(Fields &field);

  private:
    std::vector<Cell *> _cells;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                       std::map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressure(Fields &field);
    bool _temperature_boundary_type ;   ////****************NOT REQUIRED ANYMORE
    virtual void apply_temperature(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};
