#include "boundary_conditions.h"
#include <vector>
#include <cmath> // For isnan/isinf
#include <iostream>

void BoundaryConditions::apply_all(std::vector<std::vector<NodeData>>& grid, const Lattice& lattice, Real lid_velocity) {
    int nx = grid.size();
    int ny = grid[0].size();
    const auto& c = lattice.get_c();
    int Q = lattice.get_Q();

    // --- Apply Bounce-Back on Walls ---
    // Iterate through all nodes. If a node is !is_fluid, apply bounce back
    // This requires knowing which populations point *into* the fluid from the wall.
    // A simpler approach often used is to do bounce-back *during* streaming:
    // If f[opp[k]] streams from a fluid node (i,j) to a wall node (iw, jw),
    // then set f[k] at (i,j) for the *next* step to the value f[opp[k]] had *before* streaming.
    // However, the current structure applies BC *after* streaming.

    // Let's stick to post-streaming bounce-back for now.
    // For each wall node, look at its fluid neighbors.
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (!grid[i][j].is_fluid) {
                // This is a wall node. Find fluid neighbors and apply bounce-back.
                // The populations f_new[k] at the wall node came from fluid neighbors f[k] at (i-cx, j-cy).
                // We need to reflect these back.
                for (int k = 0; k < Q; ++k) {
                    int prev_i = i - static_cast<int>(c[k].x);
                    int prev_j = j - static_cast<int>(c[k].y);

                    // Check if the source node was within bounds and was fluid
                    if (prev_i >= 0 && prev_i < nx && prev_j >= 0 && prev_j < ny && grid[prev_i][prev_j].is_fluid) {
                        // The population f_new[k] at wall (i,j) came from fluid (prev_i, prev_j).
                        // Bounce it back: The population f[opp[k]] at the fluid node (prev_i, prev_j)
                        // should receive the value that just arrived at the wall.
                        grid[prev_i][prev_j].f[lattice.opposite(k)] = grid[i][j].f_new[k];
                    }
                }
            }
        }
    }


    // --- Apply Lid Velocity (Top Wall, excluding corners) ---
    Vector2D lid_u = {lid_velocity, 0.0};
    int j_top = ny - 1;
    for (int i = 1; i < nx - 1; ++i) {
        if (!grid[i][j_top].is_fluid) { // Should be wall nodes
             // We need to reconstruct the unknown populations pointing into the fluid
             // using the known velocity. Use Zou-He style non-equilibrium bounce-back.
             Real rho_neighbor = grid[i][j_top-1].rho; // Use rho from adjacent fluid node

             // Calculate missing populations (2, 6, 7 for top wall)
             // f_opp = f_k + correction
             // Correction term based on density and target velocity difference from bounce-back
             Real w2 = lattice.get_w()[2]; Real w6 = lattice.get_w()[6];
             Real w5 = lattice.get_w()[5]; // Define w5
             Real cs2 = lattice.get_cs2();

             // Populations streaming *from* the wall node (i, j_top) after bounce-back would be:
             // f_4_bb = grid[i][j_top].f_new[2]; // f_new[k] holds value streamed *to* wall node
             // f_8_bb = grid[i][j_top].f_new[6];
             // f_5_bb = grid[i][j_top].f_new[7];

             // Apply Zou-He correction based on target velocity (lid_u) and neighbor density
             // f_i = f_eq(rho_neighbor, lid_u)_i + (f_i_bb - f_eq(rho_neighbor, 0)_i) ??? No, simpler form:
             // Calculate rho at wall based on known populations + guess for unknowns assuming u=0
             // Then calculate unknowns to match target u.

             // Simpler non-equilibrium bounce back:
             // f_i = f_opp_streamed_in + 2 * w_i * rho_neighbor * (c_i . u_wall) / cs^2
             grid[i][j_top-1].f[2] = grid[i][j_top].f_new[4] + 2.0 * w2 * rho_neighbor * dot(c[2], lid_u) / cs2;
             grid[i][j_top-1].f[5] = grid[i][j_top].f_new[7] + 2.0 * w5 * rho_neighbor * dot(c[5], lid_u) / cs2;
             grid[i][j_top-1].f[6] = grid[i][j_top].f_new[8] + 2.0 * w6 * rho_neighbor * dot(c[6], lid_u) / cs2;

             // Ensure positivity (crude clamp)
             if(grid[i][j_top-1].f[2] < 0) grid[i][j_top-1].f[2] = REAL_EPSILON;
             if(grid[i][j_top-1].f[5] < 0) grid[i][j_top-1].f[5] = REAL_EPSILON;
             if(grid[i][j_top-1].f[6] < 0) grid[i][j_top-1].f[6] = REAL_EPSILON;
        }
    }
}


// --- Individual BC functions (Example - might not be used directly by apply_all above) ---

void BoundaryConditions::bounce_back(NodeData& wall_node, const Lattice& lattice, const std::vector<NodeData*>& fluid_neighbors) {
    // This function is harder to use correctly in the post-streaming application style.
    // The logic is better integrated into apply_all or done during streaming.
    // For completeness, if called *before* streaming on a wall node:
    // for (int k = 0; k < lattice.get_Q(); ++k) {
    //     wall_node.f_new[lattice.opposite(k)] = wall_node.f[k];
    // }
    (void)wall_node; // Avoid unused parameter warning
    (void)lattice;
    (void)fluid_neighbors;
}

void BoundaryConditions::fixed_velocity(NodeData& boundary_node, const Lattice& lattice, const Vector2D& target_u, Real target_rho) {
     // Zou-He method (simplified example for a specific boundary orientation)
     // Assumes boundary_node is adjacent to the fluid domain.
     // Needs information about which populations are unknown (coming from outside).
     // This is complex to generalize here. The implementation in apply_all is specific to the lid.
     (void)boundary_node; // Avoid unused parameter warning
     (void)lattice;
     (void)target_u;
     (void)target_rho;
     std::cerr << "Warning: BoundaryConditions::fixed_velocity is not fully implemented." << std::endl;
}