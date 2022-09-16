#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
namespace LidDrivenCavity {
const int moving_wall_id = 10;
const int fixed_wall_id = 6;
const double wall_velocity = 1.0;
} // namespace LidDrivenCavity

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
} // namespace border

enum class cell_type {
    FLUID,
    INFLOW,
    OUTFLOW,
    COLD_FIXED_WALL,
    HOT_FIXED_WALL,
    ADIABATIC_FIXED_WALL,
    FIXED_WALL,
    MOVING_WALL,
    HIGHER_POTENTIAL_WALL,
    LOWER_POTENTIAL_WALL,
    COUPLED_BOUNDARY,
    DEFAULT
};
