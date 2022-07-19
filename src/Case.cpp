#include "Case.hpp"
#include "Communication.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <iomanip>
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

Case::Case(std::string file_name, int argn, char **args, int rank, int size) {

    _rank = rank;
    _size = size;
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu;          /* viscosity   */
    double UI = 0.0;    /* velocity x-direction */
    double VI = 0.0;    /* velocity y-direction */
    double PI;          /* pressure */
    double GX;          /* gravitation x-direction */
    double GY;          /* gravitation y-direction */
    double xlength;     /* length of the domain x-dir.*/
    double ylength;     /* length of the domain y-dir.*/
    double dt;          /* time step */
    int imax;           /* number of cells x-direction*/
    int jmax;           /* number of cells y-direction*/
    double gamma;       /* upwind differencing factor*/
    double omg;         /* relaxation factor */
    double tau;         /* safety factor for time step*/
    int itermax;        /* max. number of iterations for pressure per time step */
    double eps;         /* accuracy bound for pressure*/
    double UIN;         /* X- Inlet Velocity*/
    double VIN;         /* Y- Inlet velocity*/
    double P_out = 0.0; /* Outlet_Pressure */

    /* WALL CLUSTERS  */
    int num_of_walls;        /* Number of walls   */
    double wall_temp_3;      /*Cold wall temperature*/
    double wall_temp_4;      /*Hot wall temperature*/
    double wall_temp_5 = -1; /* Wall temperature -1 for Adiabatic Wall  */

    /* ENERGY VARIABLES*/
    double TI;    /* Initial temperature */
    double beta;  /* Thermal Expansion Coefficient */
    double alpha; /* Thermal diffusivity */

    /*ELECTRO-MAGNETIC VARIABLES*/
    double k;    /* Electric Conductivity */
    double rho;  /* Density*/
    double phi1; /* Potential at Electrode 1*/
    double phi2; /* Potential at Electrode 2*/
    double Bz;   /* Magnetic Flux Density Perpendicular to Simulation Plane*/

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
                if (var == "Pout") file >> P_out;
                if (var == "num_of_walls") file >> num_of_walls;
                if (var == "TI") file >> TI;
                if (var == "energy_eq") {
                    std::string temp;
                    file >> temp;
                    if (temp == "on") _energy_eq = true;
                }
                if (var == "beta") file >> beta;
                if (var == "alpha") file >> alpha;
                if (var == "wall_temp_3") file >> wall_temp_3;
                if (var == "wall_temp_4") file >> wall_temp_4;
                if (var == "em_eq") {
                    std::string temp;
                    file >> temp;
                    if (temp == "on") _em_eq = true;
                }
                if (var == "v1") file >> phi1;
                if (var == "v2") file >> phi2;
                if (var == "k") file >> k;
                if (var == "rho") file >> rho;
                if (var == "Bz") file >> Bz;
                if (var == "couple") {
                    std::string temp;
                    file >> temp;
                    if (temp == "on") _couple = true;
                }
                if (var == "x_offset") file >> _x_offset;
                if (var == "y_offset") file >> _y_offset;
                if (var == "x") {
                    file >> _iproc;
                    if (_iproc < 1) {
                        std::cout << "Number of domain decomposition in x-direction cannot be less than 1.\n";
                        exit(0);
                    }
                }
                if (var == "y") {
                    file >> _jproc;
                    if (_jproc < 1) {
                        std::cout << "Number of domain decomposition in y-direction cannot be less than 1.\n";
                        exit(0);
                    }
                }
            }
        }
    }
    file.close();

    if ((_iproc * _jproc) != _size) {
        if (_rank == 0) {
            std::cout << "_iproc*_jproc != _size\nThe total number of processes must be the product of number of "
                         "processes in x and y directions. "
                      << "Please input correct value. Exiting!!!\n";
            Communication::finalize();
            exit(0);
        } else {
            Communication::finalize();
            exit(0);
        }
    }
    std::map<int, double> wall_vel;
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }

#ifdef PRECICE
    _config_file_name = args[2];
#endif

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / static_cast<double>(imax);
    domain.dy = ylength / static_cast<double>(jmax);
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain, _iproc, _jproc, _rank, _size);

    if (!_energy_eq && !_em_eq) {
        _field = Fields(_grid, nu, dt, tau, UI, VI, PI, GX, GY);
    } else if (_energy_eq) {
        _field = Fields(_grid, nu, dt, tau, alpha, beta, UI, VI, PI, TI, GX, GY);
    } else {
        _field = Fields(nu, dt, tau, k, rho, Bz, UI, VI, PI, GX, GY, _grid);
    }

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;

    std::map<int, double> temp1 = {{3, wall_temp_3}};
    std::map<int, double> temp2 = {{4, wall_temp_4}};
    std::map<int, double> temp3 = {{5, wall_temp_5}};

    std::map<int, double> temp4 = {{11, phi1}};
    std::map<int, double> temp5 = {{12, phi2}};

    // Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
    }
    if (not _grid.fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
    }
    if (not _grid.cold_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.cold_fixed_wall_cells(), temp1));
    }
    if (not _grid.hot_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.hot_fixed_wall_cells(), temp2));
    }
    if (not _grid.adiabatic_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.adiabatic_fixed_wall_cells(), temp3));
    }
    if (not _grid.inflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.inflow_cells(), UIN, VIN, TI));
    }
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutflowBoundary>(_grid.outflow_cells(), P_out));
    }
    // Constructing Potential Boundaries
    if (not _grid.higher_potential_cells().empty()) {
        _potential_boundaries.push_back(std::make_unique<PotentialBoundary>(_grid.higher_potential_cells(), temp4));
    }
    if (not _grid.lower_potential_cells().empty()) {
        _potential_boundaries.push_back(std::make_unique<PotentialBoundary>(_grid.lower_potential_cells(), temp5));
    }

// Constructing Coupled Boundaries
#ifdef PRECICE
    if (not _grid.coupled_cells().empty()) {
        _coupled_boundaries.push_back(std::make_unique<CoupledBoundary>(_grid.coupled_cells()));
        _TI = TI;
    }
#endif
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

    // Log File
    std::ofstream output_file;
    std::string outputname = _dict_name + '/' + _case_name + ".log";

    if (_rank == 0) {
        output_file.open(outputname);
        printIntro(output_file);
    }

    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    double output_counter = 0.0;
    uint8_t counter = 0;    // Counter for printing values on the console
    int number_fluid_cells; // Number of Fluid Cells in Domain

    output_vtk(timestep++); // Writing intial data

    if (!_energy_eq && !_em_eq) {

        if (_rank == 0) std::cout << "ENERGY EQUATION OFF" << std::endl;

#ifdef PRECICE

        const std::string solver_name("FluidSolver");
        const std::string mesh_name("FluidMesh");

        // constructing precice object
        precice::SolverInterface precice(solver_name, _config_file_name, _rank, _size);

        int dim = precice.getDimensions();
        int meshID = precice.getMeshID("FluidMesh");
        int vertexSize = _grid.coupled_cells().size(); // number of vertices at wet surface

        // assigning coords of coupling vertices
        std::vector<double> coords(vertexSize * dim);

        int count = 0;
        for (const auto &elem : _grid.coupled_cells()) {
            int i = elem->i();
            int j = elem->j() - 1;

            coords[count] = i * _grid.dx() + _x_offset;
            coords[count + 1] = j * _grid.dy() + _y_offset;
            count += 2;
        }

        // determine coordinates
        std::vector<int> vertexIDs(vertexSize);
        precice.setMeshVertices(meshID, vertexSize, coords.data(), vertexIDs.data());

        int U1_ID = precice.getDataID("X_Velocity1", meshID);
        int V1_ID = precice.getDataID("Y_Velocity1", meshID);
        int P1_ID = precice.getDataID("Pressure1", meshID);
        int F1_ID = precice.getDataID("X_Flux1", meshID);
        int G1_ID = precice.getDataID("Y_Flux1", meshID);
        int U2_ID = precice.getDataID("X_Velocity2", meshID);
        int V2_ID = precice.getDataID("Y_Velocity2", meshID);
        int P2_ID = precice.getDataID("Pressure2", meshID);
        int F2_ID = precice.getDataID("X_Flux2", meshID);
        int G2_ID = precice.getDataID("Y_Flux2", meshID);

        std::vector<double> U1(vertexSize);
        std::vector<double> V1(vertexSize);
        std::vector<double> P1(vertexSize);
        std::vector<double> F1(vertexSize);
        std::vector<double> G1(vertexSize);
        std::vector<double> U2(vertexSize);
        std::vector<double> V2(vertexSize);
        std::vector<double> P2(vertexSize);
        std::vector<double> F2(vertexSize);
        std::vector<double> G2(vertexSize);

        // initializing precice
        double precice_dt = precice.initialize();
        dt = std::min(dt, precice_dt);
#endif

#ifdef PRECICE
        while (precice.isCouplingOngoing()) {

            if (precice.isReadDataAvailable()) {
                precice.readBlockScalarData(U1_ID, vertexSize, vertexIDs.data(), U1.data());
                precice.readBlockScalarData(V1_ID, vertexSize, vertexIDs.data(), V1.data());
                precice.readBlockScalarData(P1_ID, vertexSize, vertexIDs.data(), P1.data());
                precice.readBlockScalarData(F1_ID, vertexSize, vertexIDs.data(), F1.data());
                precice.readBlockScalarData(G1_ID, vertexSize, vertexIDs.data(), G1.data());
            }
#else
        while (t < _t_end) {
#endif
            // Apply BCs
            for (auto &i : _boundaries) {
                i->apply(_field);
            }
#ifdef PRECICE
            for (auto &i : _coupled_boundaries) {
                i->apply_dirichlet_velocity(_field, U1, V1);
            }
#endif

            // Calculate Fluxes
            _field.calculate_fluxes(_grid);
            Communication::communicate(_field.f_matrix(), _grid.domain(), _rank);
            Communication::communicate(_field.g_matrix(), _grid.domain(), _rank);

#ifdef PRECICE
            for (auto &i : _coupled_boundaries) {
                i->apply_dirichlet_flux(_field, F1, G1);
            }
#endif

            //  Calculate RHS of PPE
            _field.calculate_rs(_grid);

            // Perform SOR Iterations
            int it = 0;
            double res = 1000.;
            while (it <= _max_iter && res >= _tolerance) {
                for (auto &i : _boundaries) {
                    i->apply_pressure(_field);
                }

#ifdef PRECICE
                for (auto &i : _coupled_boundaries) {
                    i->apply_dirichlet_pressure(_field, P1);
                }
#endif

                res = _pressure_solver->solve(_field, _grid, _boundaries); // Local sum
                res = Communication::reduce_sum(res);                      // Sum reduction over all domains
                number_fluid_cells = _grid.fluid_cells().size();
                number_fluid_cells =
                    Communication::reduce_sum(number_fluid_cells); // Sum of fluid cells over all domains
                res = std::sqrt(res / number_fluid_cells);         // Final residual
                Communication::communicate(_field.p_matrix(), _grid.domain(), _rank);
                it++;
            }

            // Calculate Velocities U and V
            _field.calculate_velocities(_grid);
            // Exchange velocities
            Communication::communicate(_field.u_matrix(), _grid.domain(), _rank);
            Communication::communicate(_field.v_matrix(), _grid.domain(), _rank);

            // Generating VTK files
            output_counter += dt;
            if (output_counter >= _output_freq) {
                output_vtk(timestep++);
                output_counter = 0;
                if (_rank == 0) {
                    std::cout << "\n[" << static_cast<int>((t / _t_end) * 100) << "%"
                              << " completed] Writing Data at t=" << t << "s\n";
                }
            }

            // Writing simulation data in a log file
            if (_rank == 0) {
                output_file << std::left << "Simulation Time[s] = " << std::setw(7) << t
                            << "\tTime Step[s] = " << std::setw(7) << dt << "\tSOR Iterations = " << std::setw(3) << it
                            << "\tSOR Residual = " << std::setw(7) << res << "\n";
            }

            // Printing info and checking for errors once in 5 runs of the loop
            if (counter == 10) {
                counter = 0;
                if (_rank == 0) {
                    std::cout << std::left << "Simulation Time[s] = " << std::setw(7) << t
                              << "\tTime Step[s] = " << std::setw(7) << dt << "\tSOR Iterations = " << std::setw(3)
                              << it << "\tSOR Residual = " << std::setw(7) << res << "\n";
                }

                // Check for unphysical behaviour
                if (check_err(_field, _grid.imax(), _grid.jmax())) exit(0);
            }
            counter++;

#ifdef PRECICE
            if (precice.isWriteDataRequired(dt)) {
                _field.get_border_U(1, U2);
                _field.get_border_V(1, V2);
                _field.get_border_P(1, P2);
                _field.get_border_F(1, F2);
                _field.get_border_G(1, G2);
                precice.writeBlockScalarData(U2_ID, vertexSize, vertexIDs.data(), U2.data());
                precice.writeBlockScalarData(V2_ID, vertexSize, vertexIDs.data(), V2.data());
                precice.writeBlockScalarData(P2_ID, vertexSize, vertexIDs.data(), P2.data());
                precice.writeBlockScalarData(F2_ID, vertexSize, vertexIDs.data(), F2.data());
                precice.writeBlockScalarData(G2_ID, vertexSize, vertexIDs.data(), G2.data());
            }

            precice_dt = precice.advance(dt);
#endif
            // Updating current time
            t = t + dt;

            //  Calculate Adaptive Time step
            dt = _field.calculate_dt(_grid);
            dt = Communication::reduce_min(dt);
#ifdef PRECICE
            dt = std::min(dt, precice_dt);
#endif
        }

#ifdef PRECICE
        // Terminating precice
        precice.finalize();
#endif

    } else if (!_em_eq && _energy_eq) {
        if (_rank == 0) {
            std::cout << "ENERGY EQN ON" << std::endl;
        }

#ifdef PRECICE

        const std::string solver_name("EnergySolver");
        const std::string mesh_name("EnergyMesh");

        // constructing precice object
        precice::SolverInterface precice(solver_name, _config_file_name, _rank, _size);

        int dim = precice.getDimensions();
        int meshID = precice.getMeshID("EnergyMesh");
        int vertexSize = _grid.coupled_cells().size(); // number of vertices at wet surface

        // assigning coords of coupling vertices
        std::vector<double> coords(vertexSize * dim);

        int count = 0;
        for (const auto &elem : _grid.coupled_cells()) {
            int i = elem->i();
            int j = elem->j() - 1;

            coords[count] = i * _grid.dx() + _x_offset;
            coords[count + 1] = j * _grid.dy() + _y_offset;
            count += 2;
        }

        // determine coordinates
        std::vector<int> vertexIDs(vertexSize);
        precice.setMeshVertices(meshID, vertexSize, coords.data(), vertexIDs.data());

        int U1_ID = precice.getDataID("X_Velocity1", meshID);
        int V1_ID = precice.getDataID("Y_Velocity1", meshID);
        int P1_ID = precice.getDataID("Pressure1", meshID);
        int F1_ID = precice.getDataID("X_Flux1", meshID);
        int G1_ID = precice.getDataID("Y_Flux1", meshID);
        int U2_ID = precice.getDataID("X_Velocity2", meshID);
        int V2_ID = precice.getDataID("Y_Velocity2", meshID);
        int P2_ID = precice.getDataID("Pressure2", meshID);
        int F2_ID = precice.getDataID("X_Flux2", meshID);
        int G2_ID = precice.getDataID("Y_Flux2", meshID);

        std::vector<double> U1(vertexSize);
        std::vector<double> V1(vertexSize);
        std::vector<double> P1(vertexSize);
        std::vector<double> F1(vertexSize);
        std::vector<double> G1(vertexSize);
        std::vector<double> U2(vertexSize);
        std::vector<double> V2(vertexSize);
        std::vector<double> P2(vertexSize);
        std::vector<double> F2(vertexSize);
        std::vector<double> G2(vertexSize);

        // initializing precice
        double precice_dt = precice.initialize();
        dt = std::min(dt, precice_dt);

        while (precice.isCouplingOngoing()) {

            if (precice.isReadDataAvailable()) {
                precice.readBlockScalarData(U1_ID, vertexSize, vertexIDs.data(), U1.data());
                precice.readBlockScalarData(V1_ID, vertexSize, vertexIDs.data(), V1.data());
                precice.readBlockScalarData(P1_ID, vertexSize, vertexIDs.data(), P1.data());
                precice.readBlockScalarData(F1_ID, vertexSize, vertexIDs.data(), F1.data());
                precice.readBlockScalarData(G1_ID, vertexSize, vertexIDs.data(), G1.data());
            }
#else
        while (t < _t_end) {
#endif
            // Apply BCs
            for (auto &i : _boundaries) {
                i->apply(_field);
                i->apply_temperature(_field);
            }
#ifdef PRECICE
            for (auto &i : _coupled_boundaries) {
                i->apply_dirichlet_velocity(_field, U1, V1);
                i->apply_temperature(_field, _TI);
            }
#endif

            // Calculate Temperatures
            _field.calculate_temperatures(_grid);
            Communication::communicate(_field.t_matrix(), _grid.domain(), _rank);

            // Calculate Fluxes
            _field.calculate_fluxes(_grid, 1);
            Communication::communicate(_field.f_matrix(), _grid.domain(), _rank);
            Communication::communicate(_field.g_matrix(), _grid.domain(), _rank);

#ifdef PRECICE
            for (auto &i : _coupled_boundaries) {
                i->apply_dirichlet_flux(_field, F1, G1);
            }
#endif

            // Calculate RHS of PPE
            _field.calculate_rs(_grid);

            // Perform SOR Iterations
            int it = 0;
            double res = 1000.;
            while (it <= _max_iter && res >= _tolerance) {
                for (auto &i : _boundaries) {
                    i->apply_pressure(_field);
                }

#ifdef PRECICE
                for (auto &i : _coupled_boundaries) {
                    i->apply_dirichlet_pressure(_field, P1);
                }
#endif

                res = _pressure_solver->solve(_field, _grid, _boundaries);

                // Sum reduction over all domains
                res = Communication::reduce_sum(res);

                number_fluid_cells = _grid.fluid_cells().size();
                number_fluid_cells =
                    Communication::reduce_sum(number_fluid_cells); // Sum of fluid cells over all domains
                res = std::sqrt(res / number_fluid_cells);         // Final residual

                Communication::communicate(_field.p_matrix(), _grid.domain(), _rank); // communicate pressures
                it++;
            }

            // Calculate Velocities U and V
            _field.calculate_velocities(_grid);
            Communication::communicate(_field.u_matrix(), _grid.domain(), _rank);
            Communication::communicate(_field.v_matrix(), _grid.domain(), _rank);

            // Storing the values in the VTK file
            output_counter += dt;
            if (output_counter >= _output_freq) {
                output_vtk(timestep++);
                output_counter = 0;
                if (_rank == 0) {
                    std::cout << "\n[" << static_cast<int>((t / _t_end) * 100) << "%"
                              << " completed] Writing Data at t=" << t << "s\n"
                              << "\n\n";
                }
            }

            // Writing simulation data in a log file
            if (_rank == 0) {
                output_file << std::left << "Simulation Time[s] = " << std::setw(7) << t
                            << "\tTime Step[s] = " << std::setw(7) << dt << "\tSOR Iterations = " << std::setw(3) << it
                            << "\tSOR Residual = " << std::setw(7) << res << "\n";
            }

            // Printing info and checking for errors once in 10 runs of the loop
            if (counter == 10) {
                counter = 0;
                if (_rank == 0) {
                    std::cout << std::left << "Simulation Time[s] = " << std::setw(7) << t
                              << "\tTime Step[s] = " << std::setw(7) << dt << "\tSOR Iterations = " << std::setw(3)
                              << it << "\tSOR Residual = " << std::setw(7) << res << "\n";
                }

                // Check for unphysical behaviour
                if (check_err(_field, _grid.imax(), _grid.jmax())) exit(0);
            }
            counter++;

#ifdef PRECICE
            if (precice.isWriteDataRequired(dt)) {
                _field.get_border_U(1, U2);
                _field.get_border_V(1, V2);
                _field.get_border_P(1, P2);
                _field.get_border_F(1, F2);
                _field.get_border_G(1, G2);
                precice.writeBlockScalarData(U2_ID, vertexSize, vertexIDs.data(), U2.data());
                precice.writeBlockScalarData(V2_ID, vertexSize, vertexIDs.data(), V2.data());
                precice.writeBlockScalarData(P2_ID, vertexSize, vertexIDs.data(), P2.data());
                precice.writeBlockScalarData(F2_ID, vertexSize, vertexIDs.data(), F2.data());
                precice.writeBlockScalarData(G2_ID, vertexSize, vertexIDs.data(), G2.data());
            }
            precice_dt = precice.advance(dt);
#endif
            // Updating current time
            t = t + dt;

            // Calculate Adaptive Time step
            dt = _field.calculate_dt_e(_grid);
            dt = Communication::reduce_min(dt);
#ifdef PRECICE
            dt = std::min(dt, precice_dt);
#endif
        }

#ifdef PRECICE
    // finalizing precice
    precice.finalize();
#endif  
    } else {

        if (_rank == 0) std::cout << "ELECTROMAGNETIC EQUATION ON" << std::endl;

#ifdef PRECICE
        const std::string solver_name("EM_Solver");
        const std::string mesh_name("EM_Pump-Mesh");

        // constructing precice object
        precice::SolverInterface precice(solver_name, _config_file_name, _rank, _size);

        int dim = precice.getDimensions();
        int meshID = precice.getMeshID("EM_Pump-Mesh");
        int vertexSize = _grid.coupled_cells().size(); // number of vertices at wet surface
        // determine vertexSize
        std::vector<double> coords(vertexSize * dim); // coords of coupling vertices

        int count = 0;
        for (const auto &elem : _grid.coupled_cells()) {
            int i = elem->i() - 1;
            int j = elem->j() - 1;

            coords[count] = i * _grid.dx() + _x_offset;
            coords[count + 1] = j * _grid.dy() + _y_offset;
            count += 2;
        }

        // determine coordinates
        std::vector<int> vertexIDs(vertexSize);
        precice.setMeshVertices(meshID, vertexSize, coords.data(), vertexIDs.data());

        int U1_ID = precice.getDataID("X_Velocity1", meshID);
        int V1_ID = precice.getDataID("Y_Velocity1", meshID);
        int P1_ID = precice.getDataID("Pressure1", meshID);
        int F1_ID = precice.getDataID("X_Flux1", meshID);
        int G1_ID = precice.getDataID("Y_Flux1", meshID);
        int U2_ID = precice.getDataID("X_Velocity2", meshID);
        int V2_ID = precice.getDataID("Y_Velocity2", meshID);
        int P2_ID = precice.getDataID("Pressure2", meshID);
        int F2_ID = precice.getDataID("X_Flux2", meshID);
        int G2_ID = precice.getDataID("Y_Flux2", meshID);

        std::vector<double> U1(vertexSize);
        std::vector<double> V1(vertexSize);
        std::vector<double> P1(vertexSize);
        std::vector<double> F1(vertexSize);
        std::vector<double> G1(vertexSize);
        std::vector<double> U2(vertexSize);
        std::vector<double> V2(vertexSize);
        std::vector<double> P2(vertexSize);
        std::vector<double> F2(vertexSize);
        std::vector<double> G2(vertexSize);

#endif

        // Solve for Potential
        double res = 1000.;
        while (res >= 1e-6) {

            for (auto &i : _potential_boundaries) {
                i->apply_potential(_field);
            }

            res = _pressure_solver->solve_potential(_field, _grid, _potential_boundaries); // Local sum
            res = Communication::reduce_sum(res); // Sum reduction over all domains
            number_fluid_cells = _grid.fluid_cells().size();
            number_fluid_cells = Communication::reduce_sum(number_fluid_cells); // Sum of fluid cells over all domains
            res = std::sqrt(res / number_fluid_cells);                          // Final residual
            Communication::communicate(_field.phi_matrix(), _grid.domain(), _rank);
        }

        // Calculate Electric Fields
        _field.calculate_electric_fields(_grid);
        Communication::communicate(_field.ex_matrix(), _grid.domain(), _rank);
        Communication::communicate(_field.ey_matrix(), _grid.domain(), _rank);

        _field.calculate_em_forces(_grid);
        Communication::communicate(_field.fx_matrix(), _grid.domain(), _rank);
        Communication::communicate(_field.fy_matrix(), _grid.domain(), _rank);

#ifdef PRECICE
        // initializing precice
        double precice_dt = precice.initialize();
        dt = std::min(dt, precice_dt);

        while (precice.isCouplingOngoing()) {
            if (precice.isReadDataAvailable()) {
                precice.readBlockScalarData(U2_ID, vertexSize, vertexIDs.data(), U2.data());
                precice.readBlockScalarData(V2_ID, vertexSize, vertexIDs.data(), V2.data());
                precice.readBlockScalarData(P2_ID, vertexSize, vertexIDs.data(), P2.data());
                precice.readBlockScalarData(F2_ID, vertexSize, vertexIDs.data(), F2.data());
                precice.readBlockScalarData(G2_ID, vertexSize, vertexIDs.data(), G2.data());
            }
#else
        while (t < _t_end) {
#endif
            // Apply BCs
            for (auto &i : _boundaries) {
                i->apply(_field);
            }

#ifdef PRECICE
            for (auto &i : _coupled_boundaries) {
                i->apply_dirichlet_velocity(_field, U2, V2);
            }
#endif

            // Calculate Fluxes
            _field.calculate_fluxes(_grid, 2);
            Communication::communicate(_field.f_matrix(), _grid.domain(), _rank);
            Communication::communicate(_field.g_matrix(), _grid.domain(), _rank);

#ifdef PRECICE
            for (auto &i : _coupled_boundaries) {
                i->apply_dirichlet_flux(_field, F2, G2);
            }
#endif

            //  Calculate RHS of PPE
            _field.calculate_rs(_grid);

            // Perform SOR Iterations
            int it = 0;
            res = 1000.;
            while (it <= _max_iter && res >= _tolerance) {
                for (auto &i : _boundaries) {
                    i->apply_pressure(_field);
                }

#ifdef PRECICE
                for (auto &i : _coupled_boundaries) {
                    i->apply_dirichlet_pressure(_field, P2);
                }
#endif

                res = _pressure_solver->solve(_field, _grid, _boundaries); // Local sum
                res = Communication::reduce_sum(res);                      // Sum reduction over all domains
                number_fluid_cells = _grid.fluid_cells().size();
                number_fluid_cells =
                    Communication::reduce_sum(number_fluid_cells); // Sum of fluid cells over all domains
                res = std::sqrt(res / number_fluid_cells);         // Final residual
                Communication::communicate(_field.p_matrix(), _grid.domain(), _rank);
                it++;
            }

            // Calculate Velocities U and V
            _field.calculate_velocities(_grid);
            // Exchange velocities
            Communication::communicate(_field.u_matrix(), _grid.domain(), _rank);
            Communication::communicate(_field.v_matrix(), _grid.domain(), _rank);

            // Generating VTK files
            output_counter += dt;
            if (output_counter >= _output_freq) {
                output_vtk(timestep++);
                output_counter = 0;
                if (_rank == 0) {
                    std::cout << "\n[" << static_cast<int>((t / _t_end) * 100) << "%"
                              << " completed] Writing Data at t=" << t << "s\n";
                }
            }

            // Writing simulation data in a log file
            if (_rank == 0) {
                output_file << std::left << "Simulation Time[s] = " << std::setw(7) << t
                            << "\tTime Step[s] = " << std::setw(7) << dt << "\tSOR Iterations = " << std::setw(3) << it
                            << "\tSOR Residual = " << std::setw(7) << res << "\n";
            }

            // Printing info and checking for errors once in 5 runs of the loop
            if (counter == 10) {
                counter = 0;
                if (_rank == 0) {
                    std::cout << std::left << "Simulation Time[s] = " << std::setw(7) << t
                              << "\tTime Step[s] = " << std::setw(7) << dt << "\tSOR Iterations = " << std::setw(3)
                              << it << "\tSOR Residual = " << std::setw(7) << res << "\n";
                }

                // Check for unphysical behaviour
                if (check_err(_field, _grid.imax(), _grid.jmax())) exit(0);
            }
            counter++;

#ifdef PRECICE
            if (precice.isWriteDataRequired(dt)) {
                _field.get_border_U(_grid.imax(), U1);
                _field.get_border_U(_grid.imax(), U1);
                _field.get_border_V(_grid.imax(), V1);
                _field.get_border_P(_grid.imax(), P1);
                _field.get_border_F(_grid.imax(), F1);
                _field.get_border_G(_grid.imax(), G1);

                precice.writeBlockScalarData(U1_ID, vertexSize, vertexIDs.data(), U1.data());
                precice.writeBlockScalarData(V1_ID, vertexSize, vertexIDs.data(), V1.data());
                precice.writeBlockScalarData(P1_ID, vertexSize, vertexIDs.data(), P1.data());
                precice.writeBlockScalarData(F1_ID, vertexSize, vertexIDs.data(), F1.data());
                precice.writeBlockScalarData(G1_ID, vertexSize, vertexIDs.data(), G1.data());
            }
#endif

            // Updating current time
            t = t + dt;

#ifdef PRECICE
            precice_dt = precice.advance(dt);
#endif

            //  Calculate Adaptive Time step
            dt = _field.calculate_dt(_grid);
            dt = Communication::reduce_min(dt);

#ifdef PRECICE
            dt = std::min(dt, precice_dt);
#endif
        }

#ifdef PRECICE
        precice.finalize();
#endif
    }

    // Storing values at the last time step
    output_vtk(timestep);
    if (_rank == 0) output_file.close();
}

void Case::output_vtk(int timestep) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().imin * dx + _x_offset;
    double y = _grid.domain().jmin * dy + _y_offset;

    { y += dy; }
    { x += dx; }

    double z = 0;

    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx + _x_offset;
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

    std::vector<vtkIdType> fixed_wall_cells;
    for (int i = 1; i <= _grid.imax(); i++) {
        for (int j = 1; j <= _grid.jmax(); j++) {
            if (_grid.cell(i, j).wall_id() != 0) {
                fixed_wall_cells.push_back(i - 1 + (j - 1) * _grid.imax());
            }
        }
    }

    for (auto t = 0; t < fixed_wall_cells.size(); t++) {
        structuredGrid->BlankCell(fixed_wall_cells.at(t));
    }

    // Pressure Array
    vtkDoubleArray *Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Temperature Array
    vtkDoubleArray *Temperature = vtkDoubleArray::New();
    Temperature->SetName("temperature");
    Temperature->SetNumberOfComponents(1);

    // Potential Array
    vtkDoubleArray *Potential = vtkDoubleArray::New();
    Potential->SetName("electric_potential");
    Potential->SetNumberOfComponents(1);

    // Electric field Array
    vtkDoubleArray *EField = vtkDoubleArray::New();
    EField->SetName("electric_field");
    EField->SetNumberOfComponents(3);

    // EM Force Array
    vtkDoubleArray *EMForce = vtkDoubleArray::New();
    EMForce->SetName("em_force");
    EMForce->SetNumberOfComponents(3);

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

    if (_energy_eq == true) {
        // Print Temperature from bottom to top
        for (int j = 1; j < _grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _grid.domain().size_x + 1; i++) {
                double temperature = _field.t(i, j);
                Temperature->InsertNextTuple(&temperature);
            }
        }

        // Add Temperature to Structured Grid
        structuredGrid->GetCellData()->AddArray(Temperature);
    }

    if (_em_eq == true) {
        float e_field[3], em_force[3];
        e_field[2] = 0;
        em_force[2] = 0;
        // Print Potential, e_field and em_force array
        for (int j = 0; j < _grid.domain().size_y + 1; j++) {
            for (int i = 0; i < _grid.domain().size_x + 1; i++) {
                e_field[0] = _field.ex(i, j);
                e_field[1] = _field.ey(i, j);
                em_force[0] = _field.fx(i, j);
                em_force[1] = _field.fy(i, j);

                EField->InsertNextTuple(e_field);
                EMForce->InsertNextTuple(em_force);
            }
        }

        for (int j = 1; j < _grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _grid.domain().size_x + 1; i++) {
                double potential = _field.phi(i, j);
                Potential->InsertNextTuple(&potential);
            }
        }

        // Add potential, e_field and em_force to Structured Grid
        structuredGrid->GetCellData()->AddArray(Potential);
        structuredGrid->GetPointData()->AddArray(EField);
        structuredGrid->GetPointData()->AddArray(EMForce);
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {

    /// Indices of rank
    int I, J;

    /// Temporary variables used to exchange information between processes
    int imin, imax, jmin, jmax, size_x, size_y;

    if (_rank == 0) {
        // std::cout << domain.size_x << " in case " << domain.size_y << std::endl;
        for (int i = 1; i < _size; ++i) {
            I = i % _iproc + 1;
            J = i / _iproc + 1;
            imin = (I - 1) * (imax_domain / _iproc);
            imax = I * (imax_domain / _iproc) + 2;
            jmin = (J - 1) * (jmax_domain / _jproc);
            jmax = J * (jmax_domain / _jproc) + 2;
            size_x = imax_domain / _iproc;
            size_y = jmax_domain / _jproc;

            // Adding the extra cells when number of cells is not divisible by iproc and jproc
            if (I == _iproc) {
                imax = imax_domain + 2;
                size_x = imax - imin - 2;
            }

            if (J == _jproc) {
                jmax = jmax_domain + 2;
                size_y = jmax - jmin - 2;
            }

            std::array<int, 4> neighbours = {-1, -1, -1, -1};
            if (I > 1) {
                // Left neighbour
                neighbours[0] = i - 1;
            }

            if (J > 1) {
                // Bottom neighbour
                neighbours[3] = i - _iproc;
            }

            if (_iproc > 1 && I < (_size - 1) % _iproc + 1) {
                // Right neighbour
                neighbours[1] = i + 1;
            }

            if (_jproc > 1 && J < (_size - 1) / _iproc + 1) {
                // Top neighbour
                neighbours[2] = i + _iproc;
            }

            MPI_Send(&imin, 1, MPI_INT, i, 999, MPI_COMM_WORLD);
            MPI_Send(&imax, 1, MPI_INT, i, 998, MPI_COMM_WORLD);
            MPI_Send(&jmin, 1, MPI_INT, i, 997, MPI_COMM_WORLD);
            MPI_Send(&jmax, 1, MPI_INT, i, 996, MPI_COMM_WORLD);
            MPI_Send(&size_x, 1, MPI_INT, i, 995, MPI_COMM_WORLD);
            MPI_Send(&size_y, 1, MPI_INT, i, 994, MPI_COMM_WORLD);
            MPI_Send(neighbours.data(), 4, MPI_INT, i, 993, MPI_COMM_WORLD);
        }
        // For rank 0
        I = _rank % _iproc + 1;
        J = _rank / _iproc + 1;
        domain.imin = (I - 1) * (imax_domain / _iproc);
        domain.imax = I * (imax_domain / _iproc) + 2;
        domain.jmin = (J - 1) * (jmax_domain / _jproc);
        domain.jmax = J * (jmax_domain / _jproc) + 2;
        domain.size_x = imax_domain / _iproc;
        domain.size_y = jmax_domain / _jproc;
        domain.neighbours[0] = -1; // left
        domain.neighbours[1] = -1; // right
        if (_iproc > 1) domain.neighbours[1] = 1;
        domain.neighbours[2] = -1; // top
        if (_jproc > 1) domain.neighbours[2] = _iproc;
        domain.neighbours[3] = -1; // bottom
    } else {
        MPI_Status status;
        MPI_Recv(&domain.imin, 1, MPI_INT, 0, 999, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.imax, 1, MPI_INT, 0, 998, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.jmin, 1, MPI_INT, 0, 997, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.jmax, 1, MPI_INT, 0, 996, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.size_x, 1, MPI_INT, 0, 995, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.size_y, 1, MPI_INT, 0, 994, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.neighbours, 4, MPI_INT, 0, 993, MPI_COMM_WORLD, &status);
    }
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

void Case::printIntro(std::ofstream &output_file) {
    output_file << "Welcome to FluidChen Flow Solver!\nThis project wqas developed as a part of Computational "
                   "Fluid Dynamics Lab course.\n\n\n\n";
    output_file << "                                                     ..::::::..\n"
                << "                                                 .^!7?7^:::.:^^:\n"
                << "                                              .~?55YY!..:~:   .\n"
                << "                                            :7Y5GGP5YJ?77:...:^^^.\n"
                << "                                          :?Y5YG#BGG5J7~^::::^~~~:\n"
                << "                                        .7Y5YYYB&#P7:.\n"
                << "                                       :J55YYYYB#!\n"
                << "                                      ^Y5YYYYYYP~\n"
                << "                                     :Y5YYYYYYYJ.\n"
                << "                                     J5YYYYYYYJ?!..::::.\n"
                << "                                    ~YYYYYYYYYJ?7~^:::::\n"
                << "                                    7YYYYYYYYJJ?77^\n"
                << "                                    ?JJJYYYYJJJ???!\n"
                << "                                    ?JJJJJJJJJJJJJJ~\n"
                << "                                    !???JJJJJJYYYYYY^ \n"
                << "                                    ^????JJJYYYY555P5^\n"
                << "                                    .????JJJYY55PPPPGJ\n"
                << "                                     ^???JJY55PPPGGGBY\n"
                << "                                      !JJJYY5PPPGGGBP:\n"
                << "                                       ^?Y55PPPGGBGJ:\n"
                << "                                         :~7JJYJ?!:\n\n\n"
                <<

        "               ?G???7. 5J     JY    J5  YY .P5?Y57   :!~^^^~^  ?.   !~  ?!:^^:  Y7.   J\n"
                << "               JG      PY     Y5    YP  55 .B7  :GY .J.        J~:::?~  ?~.::   J:7~  J.\n"
                << "               JB777:  PY     YG.   P5  55 .B7  .PY .J.        J~:^:?~  ?!:::   J  ~7:J.\n"
                << "               JP      5P???! :Y5??Y5:  5Y .G5?J5?.  :7~^:^~^  J.   !~  ?~::::  J   :7Y\n\n\n\n";
}

void Case::printIntro() {
    std::cout << "\n\n\n\nWelcome to FluidChen Flow Solver!\nThis project was developed as a part of Computational "
                 "Fluid Dynamics Lab course.\n\n\n\n";
    std::cout << "                                                     ..::::::..\n"
              << "                                                 .^!7?7^:::.:^^:\n"
              << "                                              .~?55YY!..:~:   .\n"
              << "                                            :7Y5GGP5YJ?77:...:^^^.\n"
              << "                                          :?Y5YG#BGG5J7~^::::^~~~:\n"
              << "                                        .7Y5YYYB&#P7:.\n"
              << "                                       :J55YYYYB#!\n"
              << "                                      ^Y5YYYYYYP~\n"
              << "                                     :Y5YYYYYYYJ.\n"
              << "                                     J5YYYYYYYJ?!..::::.\n"
              << "                                    ~YYYYYYYYYJ?7~^:::::\n"
              << "                                    7YYYYYYYYJJ?77^\n"
              << "                                    ?JJJYYYYJJJ???!\n"
              << "                                    ?JJJJJJJJJJJJJJ~\n"
              << "                                    !???JJJJJJYYYYYY^ \n"
              << "                                    ^????JJJYYYY555P5^\n"
              << "                                    .????JJJYY55PPPPGJ\n"
              << "                                     ^???JJY55PPPGGGBY\n"
              << "                                      !JJJYY5PPPGGGBP:\n"
              << "                                       ^?Y55PPPGGBGJ:\n"
              << "                                         :~7JJYJ?!:\n\n\n"
              <<

        "               ?G???7. 5J     JY    J5  YY .P5?Y57   :!~^^^~^  ?.   !~  ?!:^^:  Y7.   J\n"
              << "               JG      PY     Y5    YP  55 .B7  :GY .J.        J~:::?~  ?~.::   J:7~  J.\n"
              << "               JB777:  PY     YG.   P5  55 .B7  .PY .J.        J~:^:?~  ?!:::   J  ~7:J.\n"
              << "               JP      5P???! :Y5??Y5:  5Y .G5?J5?.  :7~^:^~^  J.   !~  ?~::::  J   :7Y\n\n\n\n";
}
