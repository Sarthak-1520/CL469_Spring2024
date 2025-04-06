#pragma once
#include "common_types.h"
#include "lattice.h"
#include <vector>

class NodeData {
public:
    NodeData(int q_vel); // Constructor needs Q

    // Distribution functions (f_i)
    std::vector<Real> f;
    std::vector<Real> f_eq; // Store equilibrium separately
    std::vector<Real> f_new; // For streaming update

    // Macroscopic variables
    Real rho = 1.0;
    Vector2D u = {0.0, 0.0};
    Vector2D u_old = {0.0, 0.0}; // To track convergence

    bool is_fluid = true; // Flag for boundary handling

    // Calculate rho and u from f
    void compute_macroscopics(const Lattice& lattice);
    // Initialize f based on f_eq for given rho, u
    void initialize_equilibrium(const Lattice& lattice);
};