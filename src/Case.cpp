#include "Case.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <chrono>
#include <iterator>
#include <string>

#ifdef GCC_VERSION_9_OR_HIGHER
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include <fstream>
#include <ios>
#include <iostream>
#include <map>
#include <vector>

#ifdef GCC_VERSION_9_OR_HIGHER
namespace filesystem = std::filesystem;
#else
namespace filesystem = std::experimental::filesystem;
#endif

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

Case::Case(std::string file_name, int argn, char **args) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu;       /* viscosity   */
    double UI = 0.0; /* velocity x-direction */
    double VI = 0.0; /* velocity y-direction */
    double PI;       /* pressure */
    double GX;       /* gravitation x-direction */
    double GY;       /* gravitation y-direction */
    double xlength;  /* length of the domain x-dir.*/
    double ylength;  /* length of the domain y-dir.*/
    double dt;       /* time step */
    int imax;        /* number of cells x-direction*/
    int jmax;        /* number of cells y-direction*/
    double gamma;    /* uppwind differencing factor*/
    double omg;      /* relaxation factor */
    double tau;      /* safety factor for time step*/
    int itermax;     /* max. number of iterations for pressure per time step */
    double eps;      /* accuracy bound for pressure*/

    double UIN; /* X- Inlet Velocity*/
    double VIN; /* Y- Inlet velocity*/

    double POUT = 0;  /* Outflow Pressure (Default Initialized to 0)
 
       /* WALL CLUSTERS  */
    int num_of_walls; /* Number of walls   */
    double wall_temp_3 = -1;
    double wall_temp_4 = -1;
    double wall_temp_5 = -1; /* Wall temperature -1 for Adiabatic Wall  */

    /* ENERGY VARIABLES*/
    double TI;              /* Initial temperature*/
    bool energy_eq = false; /* Set to True: Energy equation enabled*/
    double beta;            /* Thermal Expansion Coefficient  */
    double alpha;           /* Thermal diffusivity   */

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;

                if (var == "geo_file") file >> _geom_name;
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                if (var == "POUT") file >> POUT;
                if (var == "num_of_walls") file >> num_of_walls;
                if (var == "wall_temp_3") file >> wall_temp_3;
                if (var == "wall_temp_4") file >> wall_temp_4;
                if (var == "wall_temp_5") file >> wall_temp_5;
                if (var == "TI") file >> TI;
                if (var == "energy_eq") {
                    std::string temp;
                    file >> temp;
                    if (temp == "on") energy_eq = true;
                }
                if (var == "beta") file >> beta;
                if (var == "alpha") file >> alpha;
            }
        }
    }
    file.close();

    std::map<int, double> wall_vel;
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / static_cast<double>(imax);
    domain.dy = ylength / static_cast<double>(jmax);
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;

    // Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
    }
    if (not _grid.fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
    }
    if (not _grid.cold_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.cold_fixed_wall_cells(), wall_temp_3));
    }
    if (not _grid.hot_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.hot_fixed_wall_cells(), wall_temp_4));
    }
    if (not _grid.adiabatic_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.adiabatic_fixed_wall_cells(), wall_temp_5));
    }
    if (not _grid.inflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.fixed_wall_cells(), UIN, VIN));
    }
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutflowBoundary>(_grid.fixed_wall_cells(),POUT));
    }
}

void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using apply() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver class
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculat the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface. No member functions should be defined in abstract classes. You need to define functions in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate() {

    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    double output_counter = 0.0;
    uint8_t ctr = 0;

    auto start = std::chrono::steady_clock::now();

    output_vtk(t); // Writing intial data

    while (t < _t_end) {

        // Apply BCs
        for (auto &i : _boundaries) {
            i->apply(_field);
        }

        // Calculate Fluxes
        _field.calculate_fluxes(_grid);

        // Calculate RHS of PPE
        _field.calculate_rs(_grid);

        // Perform SOR Iterations
        int it = 0;
        double res = 1000.;
        while (it <= _max_iter && res >= _tolerance) {
            res = _pressure_solver->solve(_field, _grid, _boundaries);
            it++;
        }

        if (it >= _max_iter) std::cout << "\nSOR Max Iteration Reached!\nSOR Residue=" << res << "\n\n";

        // Calculate Velocities U and V
        _field.calculate_velocities(_grid);

        // Storing the values in the VTK file
        output_counter += dt;
        if (output_counter >= _output_freq) {
            output_vtk(t);
            output_counter = 0;
            std::cout << "\nWriting Data at t=" << t << "s\n\n";
        }

        // Printing info and checking for errors once in 5 runs of the loop
        if (ctr == 5) {
            ctr = 0;
            std::cout << "Simulation Time=" << t << "s         Time Step=" << dt << "s\n";

            if (check_err(_field, _grid.imax(), _grid.jmax())) exit(0); // Check for unphysical behaviour
        }
        ctr++;

        // Updating current time
        t = t + dt;

        // Calculate Adaptive Time step
        dt = _field.calculate_dt(_grid);
    }

    // Storing values at the last time step
    output_vtk(t);

    std::cout << "\nSimulation Complete!\n";
    auto end = std::chrono::steady_clock::now();
    cout << "Software Runtime:" << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s\n\n";
}

void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().imin * dx;
    double y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkDoubleArray *Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Print pressure and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = imax_domain + 2;
    domain.jmax = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
}

bool Case::check_err(Fields &field, int imax, int jmax) {
    for (int i = 0; i < imax + 2; i++) {
        for (int j = 0; j < jmax + 2; j++) {
            if (std::isnan(field.u(i, j)) || std::isinf(field.u(i, j))) {
                std::cout << "\nError!!!!!!!!!!\nValue of x-velocity at " << i << "," << j << " is:" << field.u(i, j)
                          << "\nExecution terminated!\n";
                return true;
            }

            if (std::isnan(field.v(i, j)) || std::isinf(field.v(i, j))) {
                std::cout << "\nError!!!!!!!!!!\nValue of y-velocity at " << i << "," << j << " is:" << field.v(i, j)
                          << "\nExecution terminated!\n";
                return true;
            }

            if (std::isnan(field.p(i, j)) || std::isinf(field.v(i, j))) {
                std::cout << "\nError!!!!!!!!!!\nValue of pressure at " << i << "," << j << " is:" << field.p(i, j)
                          << "\nExecution terminated!\n";
                return true;
            }
        }
    }
    return false;
}
