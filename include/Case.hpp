#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"

/**
 * @brief Class to hold and orchestrate the simulation flow.
 *
 */
class Case {
  public:
    /**
     * @brief Parallel constructor for the Case.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters Creates output directory
     *
     * @param[in] Input file name
     */
    Case(std::string file_name, int argn, char **args,int rank,int size);

    /**
     * @brief Main function to simulate the flow until the end time.
     *
     * Calculates the fluxes
     * Calculates the right hand side
     * Solves pressure
     * Calculates velocities
     * Outputs the solution files
     */
    void simulate();

    /**
     * @brief Prints introductory message
     */
    static void printIntro();

    /**
     * @brief Writes introductory message to log file
     * @param [in] Output file variable
     */
    void printIntro(std::ofstream &);

  private:
    /// Plain case name without paths
    std::string _case_name;
    /// Output directiory name
    std::string _dict_name;
    /// Geometry file name
    std::string _geom_name{"NONE"};
    /// Relative input file path
    std::string _prefix;

    /// Simulation time
    double _t_end;
    /// Solution file outputting frequency
    double _output_freq;
    Fields _field;
    Grid _grid;
    Discretization _discretization;
    std::unique_ptr<PressureSolver> _pressure_solver;
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    /// Set to true to enable energy equations
    bool _energy_eq = false;

    /// Parallelizing variables
    int _rank = 0;
    int _iproc=1;
    int _jproc=1;
    int _size=1;

    /// Solver convergence tolerance
    double _tolerance;

    /// Maximum number of iterations for the solver
    int _max_iter;

    /**
     * @brief Creating file names from given input data file
     *
     * Extracts path of the case file and creates code-readable file names
     * for outputting directory and geometry file.
     *
     * @param[in] input data file
     */
    void set_file_names(std::string file_name);

    /**
     * @brief Solution file outputter
     *
     * Outputs the solution files in .vtk format. Ghost cells are excluded.
     * Pressure is cell variable while velocity is point variable while being
     * interpolated to the cell faces
     *
     * @param[in] Timestep of the solution
     */
    void output_vtk(int t);

    void build_domain(Domain &domain, int imax_domain, int jmax_domain);

    /**
     * @brief Checks for unphysical values in velocity and pressure
     *
     * Prints the unphysical values if any and their respective index
     *
     * @param[in] Field
     * @param[in] imax
     * @param[in] jmax
     */
    bool check_err(Fields &field, int imax, int jmax);
};
