#include "boundary_conditions.h"
#include <vector>
#include <cmath> // For isnan/isinf
#include <iostream>

void BoundaryConditions::apply_all(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice,
    Real characteristic_velocity,
    GeometryType geometry_type
) {
    int nx = grid.size();
    int ny = grid[0].size();
    const auto& c = lattice.get_c();
    const auto& w = lattice.get_w(); // Need weights for inlet BC
    int Q = lattice.get_Q();
    Real cs2 = lattice.get_cs2();

    // Apply different boundary conditions based on geometry type
    switch (geometry_type) {
        case GeometryType::ELLIPSE:
        case GeometryType::CHANNEL_OBSTACLE:
            // --- Apply Bounce-Back on ALL Walls (Top/Bottom/Obstacle) ---
            apply_bounce_back_walls(grid, lattice);

            // --- Apply Inlet Condition (Left Wall, i=0) --- Fixed Velocity (Zou-He style)
            apply_inlet_boundary(grid, lattice, characteristic_velocity);

            // --- Apply Outlet Condition (Right Wall, i=nx-1) --- Simple Extrapolation
            apply_outlet_boundary(grid, lattice);
            break;

        case GeometryType::LID_DRIVEN_CAVITY:
            // --- Apply Bounce-Back on walls ---
            apply_bounce_back_walls(grid, lattice);

            // --- Apply Lid Velocity (Top Wall, j=ny-1) ---
            apply_lid_boundary(grid, lattice, characteristic_velocity, nx, ny);
            break;

        case GeometryType::BACKWARD_FACING_STEP:
            // --- Apply Bounce-Back on walls ---
            apply_bounce_back_walls(grid, lattice);

            // --- Apply Inlet Condition (Left Wall, i=0, upper part) ---
            apply_partial_inlet_boundary(grid, lattice, characteristic_velocity, nx, ny);

            // --- Apply Outlet Condition (Right Wall, i=nx-1) ---
            apply_outlet_boundary(grid, lattice);
            break;

        case GeometryType::TAYLOR_GREEN_VORTEX:
            // --- Apply Periodic Boundaries ---
            apply_periodic_boundaries(grid, lattice, nx, ny);
            break;

        case GeometryType::POISEUILLE_FLOW:
            // --- Apply Bounce-Back on walls ---
            apply_bounce_back_walls(grid, lattice);

            // --- Apply Pressure Boundary Conditions ---
            apply_pressure_boundaries(grid, lattice, nx, ny);
            break;

        default:
            std::cerr << "Warning: Unknown geometry type in boundary conditions!" << std::endl;
            // Apply default bounce-back
            apply_bounce_back_walls(grid, lattice);
    }
}

// Helper method to apply bounce-back on all wall nodes
void BoundaryConditions::apply_bounce_back_walls(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice
) {
    int nx = grid.size();
    int ny = grid[0].size();
    const auto& c = lattice.get_c();
    int Q = lattice.get_Q();

    // Iterate through all nodes. If a node is !is_fluid, apply bounce back
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (!grid[i][j].is_fluid) {
                // This is a wall node
                for (int k = 0; k < Q; ++k) {
                    int prev_i = i - static_cast<int>(c[k].x);
                    int prev_j = j - static_cast<int>(c[k].y);

                    // Check if the source node was within bounds and was fluid
                    if (prev_i >= 0 && prev_i < nx && prev_j >= 0 && prev_j < ny && grid[prev_i][prev_j].is_fluid) {
                        // Bounce it back: The population f[opp[k]] at the fluid node (prev_i, prev_j)
                        // should receive the value f_new[k] that just streamed *to* the wall node (i,j).
                        grid[prev_i][prev_j].f[lattice.opposite(k)] = grid[i][j].f_new[k];
                    }
                }
            }
        }
    }
}

// Helper method to apply inlet boundary condition
void BoundaryConditions::apply_inlet_boundary(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice,
    Real inlet_velocity
) {
    int nx = grid.size();
    int ny = grid[0].size();

    Vector2D inlet_u = {inlet_velocity, 0.0};
    Real inlet_rho = 1.0;

    for (int j = 1; j < ny - 1; ++j) { // Exclude corners (already bounce-back)
        if (grid[0][j].is_fluid) { // Check if it's a fluid node
            // Classic Zou-He for left wall inlet (velocity specified)
            Real f0 = grid[0][j].f[0]; // Use post-streamed values for known directions
            Real f2 = grid[0][j].f[2];
            Real f4 = grid[0][j].f[4];
            Real f3 = grid[0][j].f[3];
            Real f6 = grid[0][j].f[6];
            Real f7 = grid[0][j].f[7];

            // Use target density
            Real rho = inlet_rho;

            // Calculate unknown populations f1, f5, f8 to match rho and u_in
            Real ux = inlet_u.x;
            Real uy = inlet_u.y; // Should be 0

            grid[0][j].f[1] = f3 + (2.0/3.0) * rho * ux;
            grid[0][j].f[5] = f7 + 0.5 * (rho * ux + rho * uy) + (1.0/6.0) * rho * ux - 0.5 * (f2 - f4);
            grid[0][j].f[8] = f6 + 0.5 * (rho * ux - rho * uy) + (1.0/6.0) * rho * ux + 0.5 * (f2 - f4);

            // Simple positivity clamp
            if(grid[0][j].f[1] < 0) grid[0][j].f[1] = REAL_EPSILON;
            if(grid[0][j].f[5] < 0) grid[0][j].f[5] = REAL_EPSILON;
            if(grid[0][j].f[8] < 0) grid[0][j].f[8] = REAL_EPSILON;
        }
    }
}

// Helper method to apply outlet boundary condition
void BoundaryConditions::apply_outlet_boundary(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice
) {
    int nx = grid.size();
    int ny = grid[0].size();

    for (int j = 1; j < ny - 1; ++j) { // Exclude corners
        if (grid[nx-1][j].is_fluid) { // Check if it's a fluid node
            // Extrapolate populations pointing OUT of domain (leftwards: 3, 6, 7)
            grid[nx-1][j].f[3] = grid[nx-2][j].f[3];
            grid[nx-1][j].f[6] = grid[nx-2][j].f[6];
            grid[nx-1][j].f[7] = grid[nx-2][j].f[7];
        }
    }
}

// Helper method to apply lid boundary condition
void BoundaryConditions::apply_lid_boundary(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice,
    Real lid_velocity,
    int nx, int ny
) {
    Vector2D lid_u = {lid_velocity, 0.0};
    Real lid_rho = 1.0;

    for (int i = 1; i < nx - 1; ++i) { // Exclude corners
        // Top wall (j = ny-1) is the moving lid
        // Apply Zou-He BC for the lid
        Real f0 = grid[i][ny-1].f[0];
        Real f1 = grid[i][ny-1].f[1];
        Real f3 = grid[i][ny-1].f[3];
        Real f2 = grid[i][ny-1].f[2];
        Real f5 = grid[i][ny-1].f[5];
        Real f6 = grid[i][ny-1].f[6];

        // Use target density
        Real rho = lid_rho;

        // Calculate unknown populations f4, f7, f8 to match rho and lid_u
        Real ux = lid_u.x;
        Real uy = lid_u.y; // Should be 0

        grid[i][ny-1].f[4] = f2 - (2.0/3.0) * rho * uy;
        grid[i][ny-1].f[7] = f5 + 0.5 * (rho * ux - rho * uy) - (1.0/6.0) * rho * uy + 0.5 * (f1 - f3);
        grid[i][ny-1].f[8] = f6 - 0.5 * (rho * ux + rho * uy) - (1.0/6.0) * rho * uy - 0.5 * (f1 - f3);

        // Simple positivity clamp
        if(grid[i][ny-1].f[4] < 0) grid[i][ny-1].f[4] = REAL_EPSILON;
        if(grid[i][ny-1].f[7] < 0) grid[i][ny-1].f[7] = REAL_EPSILON;
        if(grid[i][ny-1].f[8] < 0) grid[i][ny-1].f[8] = REAL_EPSILON;
    }
}

// Helper method to apply partial inlet boundary condition (for backward-facing step)
void BoundaryConditions::apply_partial_inlet_boundary(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice,
    Real inlet_velocity,
    int nx, int ny
) {
    Vector2D inlet_u = {inlet_velocity, 0.0};
    Real inlet_rho = 1.0;

    // Find the step height (first fluid node from bottom)
    int step_height = 0;
    for (int j = 0; j < ny; ++j) {
        if (grid[0][j].is_fluid) {
            step_height = j;
            break;
        }
    }

    // Apply inlet BC only to the fluid nodes above the step
    for (int j = step_height; j < ny - 1; ++j) { // Exclude top corner
        if (grid[0][j].is_fluid) {
            // Classic Zou-He for left wall inlet
            Real f0 = grid[0][j].f[0];
            Real f2 = grid[0][j].f[2];
            Real f4 = grid[0][j].f[4];
            Real f3 = grid[0][j].f[3];
            Real f6 = grid[0][j].f[6];
            Real f7 = grid[0][j].f[7];

            // Use target density
            Real rho = inlet_rho;

            // Calculate unknown populations
            Real ux = inlet_u.x;
            Real uy = inlet_u.y; // Should be 0

            grid[0][j].f[1] = f3 + (2.0/3.0) * rho * ux;
            grid[0][j].f[5] = f7 + 0.5 * (rho * ux + rho * uy) + (1.0/6.0) * rho * ux - 0.5 * (f2 - f4);
            grid[0][j].f[8] = f6 + 0.5 * (rho * ux - rho * uy) + (1.0/6.0) * rho * ux + 0.5 * (f2 - f4);

            // Simple positivity clamp
            if(grid[0][j].f[1] < 0) grid[0][j].f[1] = REAL_EPSILON;
            if(grid[0][j].f[5] < 0) grid[0][j].f[5] = REAL_EPSILON;
            if(grid[0][j].f[8] < 0) grid[0][j].f[8] = REAL_EPSILON;
        }
    }
}


// Implement periodic boundary conditions
void BoundaryConditions::apply_periodic_boundaries(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice,
    int nx, int ny
) {
    const auto& c = lattice.get_c();
    int Q = lattice.get_Q();

    // Create a temporary copy of the grid to read from
    auto grid_copy = grid;

    // Apply periodic boundaries in both x and y directions
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < Q; ++k) {
                // Calculate source coordinates with periodic wrapping
                int src_i = (i - static_cast<int>(c[k].x) + nx) % nx;
                int src_j = (j - static_cast<int>(c[k].y) + ny) % ny;

                // Update distribution function
                grid[i][j].f[k] = grid_copy[src_i][src_j].f_new[k];
            }
        }
    }
}

// Implement pressure boundary conditions for Poiseuille flow
void BoundaryConditions::apply_pressure_boundaries(
    std::vector<std::vector<NodeData>>& grid,
    const Lattice& lattice,
    int nx, int ny
) {
    // Set pressure (density) at inlet and outlet
    Real inlet_rho = 1.01; // Higher pressure at inlet
    Real outlet_rho = 0.99; // Lower pressure at outlet

    // Apply pressure BC at inlet (left wall, i=0)
    for (int j = 1; j < ny - 1; ++j) { // Exclude corners
        if (grid[0][j].is_fluid) {
            pressure_boundary(grid[0][j], lattice, inlet_rho, true); // true = inlet
        }
    }

    // Apply pressure BC at outlet (right wall, i=nx-1)
    for (int j = 1; j < ny - 1; ++j) { // Exclude corners
        if (grid[nx-1][j].is_fluid) {
            pressure_boundary(grid[nx-1][j], lattice, outlet_rho, false); // false = outlet
        }
    }
}

// Implement pressure boundary condition
void BoundaryConditions::pressure_boundary(
    NodeData& boundary_node,
    const Lattice& lattice,
    Real target_rho,
    bool is_inlet
) {
    const auto& c = lattice.get_c();
    int Q = lattice.get_Q();

    // Calculate current density and velocity
    Real rho_current = 0.0;
    Vector2D u_current = {0.0, 0.0};

    // Sum known populations
    for (int k = 0; k < Q; ++k) {
        rho_current += boundary_node.f[k];
        u_current.x += c[k].x * boundary_node.f[k];
        u_current.y += c[k].y * boundary_node.f[k];
    }

    // Adjust velocity to match target density
    Real rho_diff = target_rho - rho_current;

    if (is_inlet) {
        // Inlet (left wall): adjust x-velocity
        u_current.x = 1.0 - rho_current / target_rho;
    } else {
        // Outlet (right wall): adjust x-velocity
        u_current.x = rho_current / target_rho - 1.0;
    }

    // Normalize velocity
    u_current.x /= target_rho;
    u_current.y /= target_rho;

    // Apply fixed velocity BC with target density
    fixed_velocity(boundary_node, lattice, u_current, target_rho);
}

// --- Individual BC functions (Keep for reference, but apply_all handles logic) ---

void BoundaryConditions::bounce_back(
    NodeData& wall_node,
    const Lattice& lattice,
    const std::vector<NodeData*>& fluid_neighbors
) {
    // Logic moved into apply_bounce_back_walls
    (void)wall_node; // Avoid unused parameter warning
    (void)lattice;
    (void)fluid_neighbors;
}

void BoundaryConditions::fixed_velocity(
    NodeData& boundary_node,
    const Lattice& lattice,
    const Vector2D& target_u,
    Real target_rho
) {
    // This is a simplified implementation
    // In a real implementation, we would need to handle different boundary orientations

    // For now, just set the macroscopic values
    boundary_node.rho = target_rho > 0.0 ? target_rho : boundary_node.rho;
    boundary_node.u = target_u;

    // Recalculate equilibrium distributions
    for (int k = 0; k < lattice.get_Q(); ++k) {
        boundary_node.f_eq[k] = boundary_node.calculate_equilibrium(k, lattice);
    }

    // Set distributions to equilibrium (this is a simplification)
    boundary_node.f = boundary_node.f_eq;
}