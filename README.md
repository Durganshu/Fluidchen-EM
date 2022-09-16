<div align="center">
  <img width="650" height="492" src="FluidchenEMLogo.png">
</div>

Fluidchen EM is a Computational Fluid Dynamics Solver developed for the course CFD Lab; taught at TUM Informatics, Chair of Scientific Computing in Computer Science.

It is capable of running diferent types of Fluid and Electromagnetohydrdynamics problems, which include:
- Lid-driven cavity
- Plane shear flow
- The Karman Vortex Street
- Flow over a Step
- Natural Convection
- Fluid Trap
- Rayleigh-Benard Convection
- EM Pump

With the help of preCICE, the solver can also be coupled with another solver and also with itself.

The solver can run sequentially as well as in parallel (using MPI).


## Software Requirements

This code is known to work on all currently supported Ubuntu LTS versions (22.04, 20.04, 18.04).
In particular, you will need:

- A recent version of the GCC compiler. Other compilers should also work, but you may need to tweak the CMakeLists.txt file (contributions welcome). GCC 7.4, 9.3, and 11.2 are known to work. See `CMakeLists.txt` and `src/Case.cpp` for some compiler-specific code.
- CMake, to configure the build.
- The VTK library, to generate result files. libvtk7 and libvtk9 are known to work.
- OpenMPI.
- preCICE

Get the dependencies on Ubuntu:

```shell
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install build-essential cmake libvtk7-qt-dev openmpi-bin libopenmpi-dev pkg-config git
```
For preCICE, follow these steps in succession:

```shell
wget https://github.com/precice/precice/releases/download/v2.4.0/libprecice2_2.4.0_focal.deb
sudo apt install ./libprecice2_2.4.0_focal.deb
```
If this is your first interaction with preCICE, we recommend going through the documentation on the [official website.](https://precice.org/)

## Building the code

Having installed and configured all the dependencies and wanna start quickly with the solver? Follow these steps in succession:

```shell
git clone https://gitlab.lrz.de/durganshu/group-i-cfd-lab.git
cd group-i-cfd-lab
mkdir build && cd build
cmake ..
make
```

After `make` completes successfully, you will get an executable `fluidchen_em` in your `build` directory.

### Build options

Want to do more with the solver? Well, the solver can be built in three different build configurations. Apart from the legacy `DEBUG` and `RELEASE` modes, there is one more configuation type : `PRECICE`. This can be used to couple **fluidchen_em** with any other solver (or even itself - See [preCICE](https://precice.org/)). 

By default, **fluidchen_em** is installed in `RELEASE` mode. In case you want to install in `DEBUG` or `PRECICE` mode, you can execute cmake as

```shell
cmake -DCMAKE_BUILD_TYPE= <MODE> ..
```
`<MODE>` can be `DEBUG` or `PRECICE`

You can see and modify all CMake options with, e.g., `ccmake .` inside `build/` (Ubuntu package `cmake-curses-gui`).

## Running _fluidchen_em_

### For `RELEASE` and `DEBUG` configuration types

In order to run **fluidchen_em**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. Navigate to the `build/` directory and run:

```shell
./fluidchen_em ../example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `../example_cases/LidDrivenCavity/LidDrivenCavity_Output`, which holds the `.vtk` files of the solution. The `.vtk` files can be visualized using [ParaView](https://www.paraview.org/) or any similar software.

If the input file does not contain a geometry file, `fluidchen_em` will run the lid-driven cavity case by default with the given parameters (see [below](https://gitlab.lrz.de/durganshu/group-i-cfd-lab/-/blob/EMHD_precice/README.md#geometry-files)).

### For `PRECICE` configuration type

In this case, apart from the input file, the precice configuration file (*.xml) needs to be provided as well as argument 2. An example file is present in `../example_cases/EM_Pump`. 

**Make sure to check the preCICE documentation before doing any changes to the configuration file.**

An example case is available in `/example_cases/EM_Pump/`. Open two different terminals and run these commands:

```shell
shell1$ ./fluidchen_em ../example_cases/EM_Pump/EM_Pump.dat ../example_cases/EM_Pump/precice-config.xml
shell2$ ./fluidchen_em ../example_cases/EM_Pump/FluidChannel.dat ../example_cases/EM_Pump/precice-config.xml
```
## Geometry files

We have our own simple format for reading geometry, based on the [“portable graymap” format (PGM)](https://en.wikipedia.org/wiki/Netpbm#File_formats). This is a text file with two sections: a header and the data. The header comprises of one line with a “magic number” defining its type (with P2 meaning “Portable GrayMap”, i.e. multiple gray values), one line with the number of cells per row and column (including the domain borders), as well as one line with the number of “colors” (each color being a type of cell, e.g. “fluid”, “obstacle”, “inflow”, “outflow”). After this information, a “color” for each cell follows. You can find examples in any of the subdirectories in `example_cases` (except for the LidDrivenCavity) and you can also create your own PGM files in any way you like (they are plain text files anyway).

The input files need to specify the location of the geometry files. If the input file does not contain a geometry file, `fluidchen_em` will run the lid-driven cavity case with the given parameters.
