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
    /// Elapsed time
    double elapsed_t;

    /// Ramp time interval
    double ramp_dt = 0;

    Fields() = default;

    /**
     * @brief Constructor for the fields for energy and electromagnetic equations disabled
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
     * @brief Constructor for the fields for energy equations enabled and electromagnetic equations disabled
     *
     * @param[in] grid
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] thermal diffusivity
     * @param[in] thermal expansion coefficient
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     * @param[in] initial temperature
     * @param[in] x-component of gravity
     * @param[in] y-component of gravity
     *
     */
    Fields(Grid &grid, double _nu, double _dt, double _tau, double alpha, double beta, double UI, double VI, double PI,
           double TI, double GX, double GY);

    /**
     * @brief Constructor for the fields for electromagnetic equations enabled and energy equations disabled
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] electric conductivity
     * @param[in] density
     * @param[in] magnetic flux density perpendicular to simulation plane
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     * @param[in] x-component of gravity
     * @param[in] y-component of gravity
     * @param[in] grid
     */

    Fields(double _nu, double _dt, double _tau, double _k, double _rho, double _Bz, double UI, double VI, double PI,
           double GX, double GY, Grid &grid);

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
     * @param[in] equation type:
     *  0: energy and em equation off
     *  1: energy equation on
     *  2: em equation on
     */

    void calculate_fluxes(Grid &grid, int eq_type = 0);

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
     * @brief Electric field calculation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_electric_fields(Grid &grid);

    /**
     * @brief Electromagnetic force calculation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_em_forces(Grid &grid);

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

    /// potential index based access and modify
    double &phi(int i, int j);

    /// x-electric field index based access and modify
    double &ex(int i, int j);

    /// y-electric field index based access and modify
    double &ey(int i, int j);

    /// x-electromagnetic force index based access and modify
    double &fx(int i, int j);

    /// y-electromagnetic force field index based access and modify
    double &fy(int i, int j);

    /// get timestep size
    double dt() const;

    /// field matrix access and modify
    Matrix<double> &u_matrix();
    Matrix<double> &v_matrix();
    Matrix<double> &t_matrix();
    Matrix<double> &p_matrix();
    Matrix<double> &f_matrix();
    Matrix<double> &g_matrix();
    Matrix<double> &phi_matrix();
    Matrix<double> &ex_matrix();
    Matrix<double> &ey_matrix();
    Matrix<double> &fx_matrix();
    Matrix<double> &fy_matrix();

    /// get border values
    void get_border_U(int col, std::vector<double> &U);
    void get_border_V(int col, std::vector<double> &V);
    void get_border_P(int col, std::vector<double> &P);
    void get_border_F(int col, std::vector<double> &F);
    void get_border_G(int col, std::vector<double> &G);

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
    /// electric potential
    Matrix<double> _PHI;
    /// electric field in x-direction
    Matrix<double> _Ex;
    /// electric field in y-direction
    Matrix<double> _Ey;
    /// electromagnetic force in x-direction
    Matrix<double> _Fx;
    /// electromagnetic force in y-direction
    Matrix<double> _Fy;

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
    /// electric conductivity
    double _k;
    /// density
    double _rho;
    /// magnetic flux density
    double _Bz;

    int _rank;
    int _size;
};
