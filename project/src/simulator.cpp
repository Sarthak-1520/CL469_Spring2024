#include "simulator.h"
#include <iostream>
#include <vector>
#include <cmath> // For std::sqrt, std::abs
#include <iomanip> // For std::setprecision
#include "entropic_equilibrium.h"
#include "data_writer.h" // Changed from vtk_writer.h

Simulator::Simulator(int nx_in, int ny_in, Real viscosity_in, Real lid_velocity_in, int total_steps_in, int output_freq_in)
    : nx(nx_in), ny(ny_in), viscosity(viscosity_in), lid_velocity(lid_velocity_in),
      total_steps(total_steps_in), output_freq(output_freq_in),
      lattice(), // Initialize D2Q9 lattice
      grid(nx, std::vector<NodeData>(ny, NodeData(lattice.get_Q())))
{
    // Calculate relaxation time tau from viscosity
    // viscosity = cs^2 * (tau - 0.5*dt). Assuming dt=1.
    // tau = viscosity / cs^2 + 0.5
    tau = viscosity / lattice.get_cs2() + 0.5;
    beta = 1.0 / (2.0 * tau + 1.0); // dt=1 assumed

    std::cout << "--- Simulation Parameters ---" << std::endl;
    std::cout << "Grid size: " << nx << " x " << ny << std::endl;
    std::cout << "Viscosity: " << viscosity << std::endl;
    std::cout << "Lid Velocity: " << lid_velocity << std::endl;
    std::cout << "Tau: " << tau << std::endl;
    std::cout << "Beta: " << beta << std::endl;
    std::cout << "Total steps: " << total_steps << std::endl;
    std::cout << "Output frequency: " << output_freq << std::endl;
    std::cout << "---------------------------" << std::endl;

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
    Real initial_rho = 1.0;
    Vector2D initial_u = {0.0, 0.0};

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            grid[i][j].rho = initial_rho;
            grid[i][j].u = initial_u;
            grid[i][j].u_old = initial_u;
            grid[i][j].initialize_equilibrium(lattice); // Sets f = f_eq
        }
    }
    setup_boundaries();
    std::cout << "Grid initialized." << std::endl;
}

void Simulator::setup_boundaries() {
    // Mark wall nodes (all outer boundaries except top lid)
    for (int i = 0; i < nx; ++i) {
        grid[i][0].is_fluid = false;      // Bottom wall
        if (i > 0 && i < nx - 1) {        // Exclude corners for lid velocity BC
             grid[i][ny - 1].is_fluid = false; // Top wall (lid)
        } else {
             grid[i][ny - 1].is_fluid = false; // Top corners are walls
        }
    }
    for (int j = 1; j < ny - 1; ++j) { // Exclude already set corners
        grid[0][j].is_fluid = false;      // Left wall
        grid[nx - 1][j].is_fluid = false; // Right wall
    }
     // Ensure corners are marked as non-fluid
     grid[0][0].is_fluid = false;
     grid[nx-1][0].is_fluid = false;
     grid[0][ny-1].is_fluid = false;
     grid[nx-1][ny-1].is_fluid = false;
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
}

void Simulator::collision_step() {
    #pragma omp parallel for collapse(2) // Optional: Parallelize collision
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (grid[i][j].is_fluid) {
                // 1. Compute entropic equilibrium f_eq
                if (!EntropicEquilibrium::compute(grid[i][j], lattice)) {
                    // Handle failure - maybe revert to polynomial or stop?
                    // For now, it prints a warning inside compute()
                    // We might need a fallback here if compute returns false
                }
                // 2. Perform entropic collision
                EntropicCollision::collide(grid[i][j], lattice, beta);
            }
        }
    }
}

void Simulator::streaming_step() {
    const auto& c = lattice.get_c();
    int Q = lattice.get_Q();

    #pragma omp parallel for collapse(2) // Optional: Parallelize streaming prep
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
             // Prepare f_new based on post-collision f
             // This copy is needed for the pull scheme if not done carefully in-place
             grid[i][j].f_new = grid[i][j].f;
        }
    }


    #pragma omp parallel for collapse(2) // Optional: Parallelize streaming update
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
     for (int i = 0; i < nx; ++i) {
         for (int j = 0; j < ny; ++j) {
             if (!grid[i][j].is_fluid) {
                 grid[i][j].f = grid[i][j].f_new;
             }
         }
     }
     // Now apply specific BC logic (bounce-back, lid velocity)
     BoundaryConditions::apply_all(grid, lattice, lid_velocity);
}


void Simulator::update_macroscopics_and_convergence() {
    Real max_delta_u_sq = 0.0;

    #pragma omp parallel for collapse(2) reduction(max:max_delta_u_sq)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (grid[i][j].is_fluid) {
                grid[i][j].u_old = grid[i][j].u; // Store previous velocity
                grid[i][j].compute_macroscopics(lattice); // Calculate new rho, u from post-BC f

                // Calculate change in velocity magnitude squared
                Real dux = grid[i][j].u.x - grid[i][j].u_old.x;
                Real duy = grid[i][j].u.y - grid[i][j].u_old.y;
                Real delta_u_sq = dux * dux + duy * duy;
                if (delta_u_sq > max_delta_u_sq) {
                    max_delta_u_sq = delta_u_sq;
                }
            } else {
                 // Reset wall node macroscopics (optional, compute_macroscopics handles it)
                 grid[i][j].rho = 0.0;
                 grid[i][j].u = {0.0, 0.0};
                 grid[i][j].u_old = {0.0, 0.0};
            }
        }
    }

    Real max_delta_u = std::sqrt(max_delta_u_sq);
    if (convergence_out.is_open()) {
         // Get current step number (needs to be passed or stored)
         // Assuming we call this once per step, need step counter 't'
         // Let's assume 't' is available or passed. For now, placeholder:
         // convergence_out << t << " " << std::scientific << std::setprecision(10) << max_delta_u << std::endl;
         // We'll write from the main loop where 't' is known.
    }
}

void Simulator::write_output(int time_step) {
    std::string filename = output_dir + "/lbm_output_" + std::to_string(time_step) + ".dat"; // Changed extension to .dat
    std::cout << "Writing output to " << filename << "..." << std::endl;
    DataWriter::write_dat(filename, grid, nx, ny); // Changed function call

    // Write convergence data here as well
    Real max_delta_u_sq = 0.0;
     for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
             if (grid[i][j].is_fluid) {
                Real dux = grid[i][j].u.x - grid[i][j].u_old.x;
                Real duy = grid[i][j].u.y - grid[i][j].u_old.y;
                max_delta_u_sq = std::max(max_delta_u_sq, dux * dux + duy * duy);
             }
        }
     }
     Real max_delta_u = std::sqrt(max_delta_u_sq);
      if (convergence_out.is_open()) {
          convergence_out << time_step << " " << std::scientific << std::setprecision(10) << max_delta_u << std::endl;
      }

}