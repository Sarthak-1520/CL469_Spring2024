#include "simulator.h"
#include <iostream>
#include <vector>
#include <cmath> // For std::sqrt, std::abs
#include <iomanip> // For std::setprecision
#include <filesystem> // For creating directories
#include "entropic_equilibrium.h"
#include "data_writer.h" // Changed from vtk_writer.h

std::string Simulator::get_case_name() const {
    switch (geometry_type) {
        case GeometryType::ELLIPSE:
            return "ellipse_flow";
        case GeometryType::LID_DRIVEN_CAVITY:
            return "lid_driven_cavity";
        case GeometryType::CHANNEL_OBSTACLE:
            return "channel_obstacle";
        case GeometryType::BACKWARD_FACING_STEP:
            return "backward_facing_step";
        case GeometryType::TAYLOR_GREEN_VORTEX:
            return "taylor_green_vortex";
        case GeometryType::POISEUILLE_FLOW:
            return "poiseuille_flow";
        default:
            return "unknown";
    }
}

std::string Simulator::get_collision_model_name() const {
    return (collision_model == CollisionModel::ENTROPIC) ? "entropic" : "trt";
}

Simulator::Simulator(int nx_in, int ny_in, Real viscosity_in, Real characteristic_velocity_in,
                     int total_steps_in, int output_freq_in, CollisionModel model, Real magic_param_in,
                     GeometryType geometry_type_in, Real param1, Real param2, Real param3, Real param4)
    : nx(nx_in), ny(ny_in), viscosity(viscosity_in), characteristic_velocity(characteristic_velocity_in),
      total_steps(total_steps_in), output_freq(output_freq_in), collision_model(model),
      magic_param(magic_param_in), geometry_type(geometry_type_in),
      geom_param1(param1), geom_param2(param2), geom_param3(param3), geom_param4(param4),
      lattice(), // Initialize D2Q9 lattice
      grid(nx, std::vector<NodeData>(ny, NodeData(lattice.get_Q())))
{
    // Calculate relaxation time tau from viscosity
    // viscosity = cs^2 * (tau - 0.5*dt). Assuming dt=1.
    // tau = viscosity / lattice.get_cs2() + 0.5
    tau = viscosity / lattice.get_cs2() + 0.5;
    beta = 1.0 / (2.0 * tau + 1.0); // dt=1 assumed

    // For TRT model, calculate tau_minus based on the magic parameter
    if (collision_model == CollisionModel::TRT) {
        tau_minus = TRTCollision::compute_tau_minus(tau, magic_param);
    }

    // Print geometry type
    std::string geometry_name;
    switch (geometry_type) {
        case GeometryType::ELLIPSE:
            geometry_name = "Flow around elliptical obstacle";
            break;
        case GeometryType::LID_DRIVEN_CAVITY:
            geometry_name = "Lid-driven cavity flow";
            break;
        case GeometryType::CHANNEL_OBSTACLE:
            geometry_name = "Channel flow with circular obstacle";
            break;
        case GeometryType::BACKWARD_FACING_STEP:
            geometry_name = "Backward-facing step flow";
            break;
        case GeometryType::TAYLOR_GREEN_VORTEX:
            geometry_name = "Taylor-Green vortex decay";
            break;
        case GeometryType::POISEUILLE_FLOW:
            geometry_name = "Poiseuille flow in a channel";
            break;
        default:
            geometry_name = "Unknown geometry";
    }

    std::cout << "--- Simulation Parameters ---" << std::endl;
    std::cout << "Geometry: " << geometry_name << std::endl;
    std::cout << "Grid size: " << nx << " x " << ny << std::endl;
    std::cout << "Viscosity: " << viscosity << std::endl;
    std::cout << "Characteristic velocity: " << characteristic_velocity << std::endl;
    std::cout << "Tau: " << tau << std::endl;
    std::cout << "Collision model: " << (collision_model == CollisionModel::ENTROPIC ? "Entropic" : "TRT") << std::endl;

    if (collision_model == CollisionModel::ENTROPIC) {
        std::cout << "Beta: " << beta << std::endl;
    } else {
        std::cout << "Magic parameter: " << magic_param << std::endl;
        std::cout << "Tau minus: " << tau_minus << std::endl;
    }

    std::cout << "Total steps: " << total_steps << std::endl;
    std::cout << "Output frequency: " << output_freq << std::endl;
    std::cout << "---------------------------" << std::endl;

    // Set up output directory and file paths
    std::string case_name = get_case_name();
    std::string model_name = get_collision_model_name();
    std::string case_dir = "results/" + case_name + "_" + model_name;

    // Create output directory if it doesn't exist
    std::filesystem::create_directories(case_dir);

    // Set output directory and convergence file path
    output_dir = case_dir;
    convergence_file = case_dir + "/convergence.dat";

    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << "Convergence file: " << convergence_file << std::endl;

    // Open convergence file
    convergence_out.open(convergence_file);
    if (!convergence_out) {
        std::cerr << "Error: Could not open convergence file: " << convergence_file << std::endl;
    } else {
        convergence_out << "# Step MaxDeltaU" << std::endl; // Header
    }
}

void Simulator::initialize_grid() {
    std::cout << "Initializing grid..." << std::endl;

    // Set all nodes as fluid initially
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            grid[i][j].is_fluid = true;
        }
    }

    // Setup boundaries based on geometry type
    setup_boundaries();

    // Initialize all nodes with appropriate initial conditions
    Real initial_rho = 1.0;
    Vector2D initial_u = {0.0, 0.0}; // Default initial velocity

    // Set initial conditions based on geometry type
    switch (geometry_type) {
        case GeometryType::ELLIPSE:
        case GeometryType::CHANNEL_OBSTACLE:
        case GeometryType::BACKWARD_FACING_STEP:
            // Uniform inlet flow profile
            initial_u = {characteristic_velocity, 0.0};
            break;

        case GeometryType::LID_DRIVEN_CAVITY:
            // Zero initial velocity inside cavity
            initial_u = {0.0, 0.0};
            break;

        case GeometryType::TAYLOR_GREEN_VORTEX:
            // Taylor-Green vortex has a specific initial condition
            // Will be set in the loop below
            break;

        case GeometryType::POISEUILLE_FLOW:
            // Zero initial velocity for Poiseuille flow
            initial_u = {0.0, 0.0};
            break;
    }

    // Initialize all nodes (including boundary nodes)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            grid[i][j].rho = initial_rho;

            // Special case for Taylor-Green vortex
            if (geometry_type == GeometryType::TAYLOR_GREEN_VORTEX) {
                // Taylor-Green vortex initial condition
                Real x = static_cast<Real>(i) / nx;
                Real y = static_cast<Real>(j) / ny;
                grid[i][j].u.x = characteristic_velocity * std::sin(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
                grid[i][j].u.y = -characteristic_velocity * std::cos(2.0 * M_PI * x) * std::sin(2.0 * M_PI * y);
            } else {
                grid[i][j].u = initial_u;
            }

            // If node was marked as wall by setup_boundaries, reset velocity
            if (!grid[i][j].is_fluid) {
                grid[i][j].u = {0.0, 0.0};
            }

            grid[i][j].u_old = grid[i][j].u; // Initialize u_old
            grid[i][j].initialize_equilibrium(lattice); // Sets f = f_eq based on local rho, u
        }
    }

    std::cout << "Grid initialized." << std::endl;
}

void Simulator::setup_boundaries() {
    std::cout << "Setting up boundaries for " <<
        (geometry_type == GeometryType::ELLIPSE ? "elliptical flow" :
         geometry_type == GeometryType::LID_DRIVEN_CAVITY ? "lid-driven cavity" :
         geometry_type == GeometryType::CHANNEL_OBSTACLE ? "channel with obstacle" :
         geometry_type == GeometryType::BACKWARD_FACING_STEP ? "backward-facing step" :
         geometry_type == GeometryType::TAYLOR_GREEN_VORTEX ? "Taylor-Green vortex" :
         geometry_type == GeometryType::POISEUILLE_FLOW ? "Poiseuille flow" : "unknown geometry")
        << "..." << std::endl;

    // Call the appropriate setup method based on geometry type
    switch (geometry_type) {
        case GeometryType::ELLIPSE:
            setup_ellipse_flow();
            break;
        case GeometryType::LID_DRIVEN_CAVITY:
            setup_lid_driven_cavity();
            break;
        case GeometryType::CHANNEL_OBSTACLE:
            setup_channel_obstacle();
            break;
        case GeometryType::BACKWARD_FACING_STEP:
            setup_backward_facing_step();
            break;
        case GeometryType::TAYLOR_GREEN_VORTEX:
            setup_taylor_green_vortex();
            break;
        case GeometryType::POISEUILLE_FLOW:
            setup_poiseuille_flow();
            break;
        default:
            std::cerr << "Error: Unknown geometry type!" << std::endl;
    }

    std::cout << "Boundaries set." << std::endl;
}

void Simulator::run() {
    initialize_grid();

    for (int t = 0; t <= total_steps; ++t) {
        time_step(t);

        if (t % 100 == 0) { // Print progress less frequently
             std::cout << "Step: " << t << "/" << total_steps << std::endl;
        }

        // Output data
        if (t % output_freq == 0) {
            write_output(t);
        }
    }
     convergence_out.close();
     std::cout << "Simulation finished." << std::endl;
}

void Simulator::time_step(int t) {
    collision_step();
    streaming_step();
    apply_boundary_conditions_step(); // Apply BCs *after* streaming
    update_macroscopics_and_convergence(); // Calculate rho, u and check convergence

    // Write convergence data at each time step
    if (convergence_out.is_open()) {
        convergence_out << t << " "
                       << std::scientific << std::setprecision(10)
                       << current_max_delta_u << std::endl;
    }
}

void Simulator::collision_step() {
    // Removed OpenMP pragma for compatibility
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (grid[i][j].is_fluid) {
                // 1. Compute equilibrium distribution
                if (!EntropicEquilibrium::compute(grid[i][j], lattice)) {
                    // Handle failure - maybe revert to polynomial or stop?
                    // For now, it prints a warning inside compute()
                }

                // 2. Perform collision based on selected model
                if (collision_model == CollisionModel::ENTROPIC) {
                    // Entropic collision with variable relaxation parameter
                    EntropicCollision::collide(grid[i][j], lattice, beta);
                } else {
                    // Two-Relaxation-Time collision
                    TRTCollision::collide(grid[i][j], lattice, tau, tau_minus);
                }
            }
        }
    }
}

void Simulator::streaming_step() {
    const auto& c = lattice.get_c();
    int Q = lattice.get_Q();

    // Removed OpenMP pragma for compatibility
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
             // Prepare f_new based on post-collision f
             // This copy is needed for the pull scheme if not done carefully in-place
             grid[i][j].f_new = grid[i][j].f;
        }
    }


    // Removed OpenMP pragma for compatibility
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            // Pull scheme: Update f[i][j] from neighbors' f_new
            for (int k = 0; k < Q; ++k) {
                // Determine source node coordinates
                int source_i = (i - static_cast<int>(c[k].x) + nx) % nx; // Periodic wrap (implicit here)
                int source_j = (j - static_cast<int>(c[k].y) + ny) % ny; // Periodic wrap (implicit here)

                // In a non-periodic domain, check bounds:
                int src_i_np = i - static_cast<int>(c[k].x);
                int src_j_np = j - static_cast<int>(c[k].y);

                if (src_i_np >= 0 && src_i_np < nx && src_j_np >= 0 && src_j_np < ny) {
                     grid[i][j].f[k] = grid[src_i_np][src_j_np].f_new[k]; // Pull from neighbor's post-collision state
                } else {
                     // Handle out-of-bounds pull (shouldn't happen if BCs applied correctly later)
                     // For now, maybe set to zero or equilibrium? Or rely on BCs to fix it.
                     // Let's assume BCs handle the populations pointing outwards.
                     // We set f[k] at the boundary node based on BCs in the next step.
                     // For fluid nodes pulling from a boundary, this value will be overwritten by BCs anyway? No.
                     // Let's rely on the post-streaming BC application.
                     // We copy f_new to f in apply_boundary_conditions for wall nodes.
                }
            }
        }
    }
     // After this loop, grid[i][j].f contains the post-streaming values *before* BCs.
     // grid[i][j].f_new still holds the post-collision values.
}


void Simulator::apply_boundary_conditions_step() {
     // Copy f_new to f for wall nodes first, as streaming doesn't update them
     // This might be needed BEFORE specific BCs if they read f[k]
     // However, standard Zou-He and bounce-back read f_new[k] (post-stream)
     // Let's keep the apply_all call first

     // Apply specific BC logic (bounce-back, inlet, outlet)
     BoundaryConditions::apply_all(grid, lattice, characteristic_velocity, geometry_type);

    // It might be necessary to re-copy f_new to f for wall nodes *after* BCs
    // if the BCs modify f directly instead of f_new for reflection.
    // The current bounce-back modifies f of the *neighboring* fluid node.
    // Inlet/Outlet modify f of the boundary node.
     for (int i = 0; i < nx; ++i) {
         for (int j = 0; j < ny; ++j) {
             if (!grid[i][j].is_fluid) {
                 grid[i][j].f = grid[i][j].f_new; // Ensure wall nodes have post-stream state for next step's collision (if it were calculated)
             }
         }
     }
}


void Simulator::update_macroscopics_and_convergence() {
    Real max_delta_u_sq = 0.0;

    // Removed OpenMP pragma for compatibility
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (grid[i][j].is_fluid) {
                // Store previous velocity before updating
                grid[i][j].u_old = grid[i][j].u;

                // Calculate new rho, u from post-BC f
                grid[i][j].compute_macroscopics(lattice);

                // Calculate change in velocity magnitude squared
                Real dux = grid[i][j].u.x - grid[i][j].u_old.x;
                Real duy = grid[i][j].u.y - grid[i][j].u_old.y;
                Real delta_u_sq = dux * dux + duy * duy;

                // Update maximum delta_u
                if (delta_u_sq > max_delta_u_sq) {
                    max_delta_u_sq = delta_u_sq;
                }
            } else {
                // Reset wall node macroscopics
                grid[i][j].rho = 0.0;
                grid[i][j].u = {0.0, 0.0};
                grid[i][j].u_old = {0.0, 0.0};
            }
        }
    }

    // Store the max_delta_u value for use in write_output
    current_max_delta_u = std::sqrt(max_delta_u_sq);
}

void Simulator::write_output(int time_step) {
    // Write simulation data to file with case-specific name
    std::string case_name = get_case_name();
    std::string model_name = get_collision_model_name();
    std::string filename = output_dir + "/" + case_name + "_" + model_name + "_" + std::to_string(time_step) + ".dat";
    std::cout << "Writing output to " << filename << "..." << std::endl;
    DataWriter::write_dat(filename, grid, nx, ny);
}