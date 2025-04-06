#pragma once
#include "node_data.h"
#include "lattice.h"
#include <vector>

class BoundaryConditions {
public:
    // Apply all boundary conditions to the grid
    static void apply_all(std::vector<std::vector<NodeData>>& grid, const Lattice& lattice, Real lid_velocity);

    // Simple bounce-back for solid walls (assumes node.is_fluid == false)
    static void bounce_back(NodeData& wall_node, const Lattice& lattice, const std::vector<NodeData*>& fluid_neighbors);

    // Fixed velocity boundary (e.g., for the lid) - Zou/He style simplified
    static void fixed_velocity(NodeData& boundary_node, const Lattice& lattice, const Vector2D& target_u, Real target_rho = -1.0); // target_rho < 0 means keep node rho
};