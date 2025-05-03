#pragma once
#include "node_data.h"
#include "lattice.h"
#include "boundary_conditions.h"
#include "entropic_collision.h"
#include "entropic_equilibrium.h"
#include "trt_collision.h"
#include "vtk_writer.h"
#include <vector>
#include <string>
#include <fstream> // For convergence data

class Simulator {
public:
    // Collision model types
    enum class CollisionModel {
        ENTROPIC,  // Entropic LBM with variable relaxation parameter
        TRT        // Two-Relaxation-Time model
    };

    // Using GeometryType from boundary_conditions.h

    /**
     * @brief Constructor for the simulator
     *
     * @param nx Grid size in x direction
     * @param ny Grid size in y direction
     * @param viscosity Fluid viscosity
     * @param characteristic_velocity Characteristic velocity (inlet or lid)
     * @param total_steps Total number of simulation steps
     * @param output_freq Frequency of output writing
     * @param model Collision model (ENTROPIC or TRT)
     * @param magic_param Magic parameter for TRT model
     * @param geometry_type Type of geometry to simulate
     * @param param1 First geometry parameter (depends on geometry type)
     * @param param2 Second geometry parameter (depends on geometry type)
     * @param param3 Third geometry parameter (depends on geometry type)
     * @param param4 Fourth geometry parameter (depends on geometry type)
     */
    Simulator(int nx, int ny, Real viscosity, Real characteristic_velocity,
              int total_steps, int output_freq,
              CollisionModel model = CollisionModel::ENTROPIC,
              Real magic_param = 0.25,
              GeometryType geometry_type = GeometryType::ELLIPSE,
              Real param1 = 0.0, Real param2 = 0.0, Real param3 = 0.0, Real param4 = 0.0);

    void run();

private:
    // Grid parameters
    int nx, ny;
    Real viscosity;
    Real characteristic_velocity; // Inlet or lid velocity
    int total_steps;
    int output_freq;
    Real tau; // Relaxation time
    Real beta; // Relaxation parameter for collision (dt=1 assumed)

    // Collision model parameters
    CollisionModel collision_model;
    Real tau_minus; // Second relaxation time for TRT
    Real magic_param; // Magic parameter for TRT

    // Geometry parameters
    GeometryType geometry_type;
    Real geom_param1, geom_param2, geom_param3, geom_param4;

    // Core components
    Lattice lattice;
    std::vector<std::vector<NodeData>> grid; // 2D grid

    // Convergence tracking
    Real current_max_delta_u = 0.0; // Stores the current maximum velocity change

    // Output
    std::string output_dir;
    std::string convergence_file;
    std::ofstream convergence_out;

    // Get case name based on geometry type
    std::string get_case_name() const;
    // Get collision model name
    std::string get_collision_model_name() const;

    // Initialization
    void initialize_grid();
    void setup_boundaries(); // Mark boundary nodes based on geometry_type

    // Geometry-specific initialization
    void setup_ellipse_flow();
    void setup_lid_driven_cavity();
    void setup_channel_obstacle();
    void setup_backward_facing_step();
    void setup_taylor_green_vortex();
    void setup_poiseuille_flow();

    // Main loop steps
    void time_step(int t);
    void collision_step();
    void streaming_step();
    void apply_boundary_conditions_step();
    void update_macroscopics_and_convergence(); // Combine for efficiency

    // Output
    void write_output(int time_step);
};