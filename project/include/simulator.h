#pragma once
#include "node_data.h"
#include "lattice.h"
#include "boundary_conditions.h"
#include "entropic_collision.h"
#include "entropic_equilibrium.h"
#include "vtk_writer.h"
#include <vector>
#include <string>
#include <fstream> // For convergence data

class Simulator {
public:
    // Constructor takes simulation parameters
    Simulator(int nx, int ny, Real viscosity, Real inlet_velocity, int total_steps, int output_freq);

    void run();

private:
    // Parameters
    int nx, ny;
    Real viscosity;
    Real inlet_velocity; // Renamed from lid_velocity
    int total_steps;
    int output_freq;
    Real tau; // Relaxation time
    Real beta; // Relaxation parameter for collision (dt=1 assumed)

    // Core components
    Lattice lattice;
    std::vector<std::vector<NodeData>> grid; // 2D grid

    // Output
    std::string output_dir = "results";
    std::string convergence_file = "results/convergence.dat";
    std::ofstream convergence_out;

    // Initialization
    void initialize_grid();
    void setup_boundaries(Real ellipse_cx, Real ellipse_cy, Real ellipse_a, Real ellipse_b); // Pass ellipse params
    void setup_boundaries(); // Mark boundary nodes

    // Main loop steps
    void time_step(int t);
    void collision_step();
    void streaming_step();
    void apply_boundary_conditions_step();
    void update_macroscopics_and_convergence(); // Combine for efficiency

    // Output
    void write_output(int time_step);
};