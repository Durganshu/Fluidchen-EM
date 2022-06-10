#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields for energy equations disabled
     *
     * @param[in] grid
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     * @param[in] x-component of gravity
     * @param[in] y-component of gravity
     */

    Fields(Grid &grid, double _nu, double _dt, double _tau, double UI, double VI, double PI, double GX, double GY);
    /**
     * @brief Constructor for the fields for energy equations enabled
     *
     * @param[in] grid
     * @param[in] kinematic viscosity
     * @param[in] thermal diffusivity
     * @param[in] thermal expansion coefficient
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     * @param[in] initial temperature
     * @param[in] x-component of gravity
     * @param[in] y-component of gravity
     *
     */
    Fields(Grid &grid, double nu, double alpha, double beta, double dt, double tau, double UI,
           double VI, double PI, double TI, double GX, double GY);

    /**
     * @brief Calculates the temperature based on explicit discretization of energy 
     * equations
     *
     * @param[in] grid in which the fluxes are calculated
     *
     */
    void calculate_temperatures(Grid &grid);

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     *
     */

    void calculate_fluxes(Grid &grid, bool energy_eq = 0);
    

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition without energy equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    double calculate_dt(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition with energy equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    double calculate_dt_e(Grid &grid);

    /// x-velocity index based access and modify
    double &u(int i, int j);

    /// y-velocity index based access and modify
    double &v(int i, int j);

    /// pressure index based access and modify
    double &p(int i, int j);

    /// temerature index based access and modify
    double &t(int i, int j);

    /// RHS index based access and modify
    double &rs(int i, int j);

    /// x-momentum flux index based access and modify
    double &f(int i, int j);

    /// y-momentum flux index based access and modify
    double &g(int i, int j);

    /// get timestep size
    double dt() const;

    /// pressure matrix access and modify
    Matrix<double> &p_matrix();

  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// y-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;
    /// temerature matrix
    Matrix<double> _T;
    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;

    /// kinematic viscosity
    double _nu;
    /// thermal diffusivity
    double _alpha;
    /// thermal expansion coefficient
    double _beta;
    /// gravitional accelearation in x direction
    double _gx{0.0};
    /// gravitional accelearation in y direction
    double _gy{0.0};
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;

    int _rank;
    int _size;
};
