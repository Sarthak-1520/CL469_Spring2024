#include "boundary_conditions.h"
#include <vector>
#include <cmath> // For isnan/isinf
#include <iostream>

void BoundaryConditions::apply_all(std::vector<std::vector<NodeData>>& grid, const Lattice& lattice, Real inlet_velocity) {
    int nx = grid.size();
    int ny = grid[0].size();
    const auto& c = lattice.get_c();
    const auto& w = lattice.get_w(); // Need weights for inlet BC
    int Q = lattice.get_Q();
    Real cs2 = lattice.get_cs2();

    // --- Apply Bounce-Back on ALL Walls (Top/Bottom/Ellipse) ---
    // Iterate through all nodes. If a node is !is_fluid, apply bounce back
    // to the populations in the neighboring fluid nodes that just streamed into the wall.
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (!grid[i][j].is_fluid) {
                // This is a wall node (ellipse or top/bottom).
                // Populations f_new[k] at this wall node came from fluid neighbors f[k] at (i-cx, j-cy).
                // Reflect these back to the fluid neighbors.
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


    // --- Apply Inlet Condition (Left Wall, i=0) --- Fixed Velocity (Zou-He style)
    Vector2D inlet_u = {inlet_velocity, 0.0};
    Real inlet_rho = 1.0;
    for (int j = 1; j < ny - 1; ++j) { // Exclude corners (already bounce-back)
        if (grid[0][j].is_fluid) { // Check if it's a fluid node (it should be)
            // Unknown populations: f1, f5, f8 (pointing rightwards)
            // Known populations (streamed from i=1): f3, f6, f7
            // Known populations (streamed from i=0): f0, f2, f4 (these are updated by collision/streaming before BC)

            // Use Zou-He method to find f1, f5, f8 based on inlet_rho and inlet_u
            // First, calculate rho based on known populations (assuming f1,f5,f8 are equilibrium for rho=1,u=0 initially? No, use target rho)
            // rho = f0 + f2 + f4 + 2*(f3 + f6 + f7) ?? No, that's for BB

            // Classic Zou-He for left wall inlet (velocity specified):
            Real f0 = grid[0][j].f[0]; // Use post-streamed values for known directions
            Real f2 = grid[0][j].f[2];
            Real f4 = grid[0][j].f[4];
            Real f3 = grid[0][j].f[3];
            Real f6 = grid[0][j].f[6];
            Real f7 = grid[0][j].f[7];

            // Calculate rho first (this differs slightly between Zou-He variants)
            // A common way: Assume rho = sum(f_known) / (1 - ux)
            // Let's use the target rho directly (inlet_rho = 1.0)
            Real rho = inlet_rho; 

            // Calculate unknown populations f1, f5, f8 to match rho and u_in
            Real ux = inlet_u.x;
            Real uy = inlet_u.y; // Should be 0

            grid[0][j].f[1] = f3 + (2.0/3.0) * rho * ux;
            grid[0][j].f[5] = f7 + 0.5 * (rho * ux + rho * uy) + (1.0/6.0) * rho * ux - 0.5 * (f2 - f4);
            grid[0][j].f[8] = f6 + 0.5 * (rho * ux - rho * uy) + (1.0/6.0) * rho * ux + 0.5 * (f2 - f4);

            // Simple positivity clamp (crude, may indicate instability if needed often)
            if(grid[0][j].f[1]<0) grid[0][j].f[1] = REAL_EPSILON;
            if(grid[0][j].f[5]<0) grid[0][j].f[5] = REAL_EPSILON;
            if(grid[0][j].f[8]<0) grid[0][j].f[8] = REAL_EPSILON;
        }
    }

    // --- Apply Outlet Condition (Right Wall, i=nx-1) --- Simple Extrapolation
    for (int j = 1; j < ny - 1; ++j) { // Exclude corners
         if (grid[nx-1][j].is_fluid) { // Check if it's a fluid node (it should be)
             // Extrapolate populations pointing OUT of domain (leftwards: 3, 6, 7)
             // from the previous node (nx-2)
             grid[nx - 1][j].f[3] = grid[nx - 2][j].f[3];
             grid[nx - 1][j].f[6] = grid[nx - 2][j].f[6];
             grid[nx - 1][j].f[7] = grid[nx - 2][j].f[7];

            // Optionally: recalculate rho/u at outlet based on all populations?
            // For simple extrapolation, often just applying it to the distributions is done.
         }
    }
}


// --- Individual BC functions (Keep for reference, but apply_all handles logic) ---

void BoundaryConditions::bounce_back(NodeData& wall_node, const Lattice& lattice, const std::vector<NodeData*>& fluid_neighbors) {
    // Logic moved into apply_all for post-streaming bounce-back
    (void)wall_node; // Avoid unused parameter warning
    (void)lattice;
    (void)fluid_neighbors;
}

void BoundaryConditions::fixed_velocity(NodeData& boundary_node, const Lattice& lattice, const Vector2D& target_u, Real target_rho) {
     // Logic implemented specifically for inlet in apply_all
     (void)boundary_node; // Avoid unused parameter warning
     (void)lattice;
     (void)target_u;
     (void)target_rho;
     // std::cerr << "Warning: BoundaryConditions::fixed_velocity is not fully implemented as generic function." << std::endl;
}