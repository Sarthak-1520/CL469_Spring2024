#include "node_data.h"
#include "entropic_equilibrium.h" // Needed for initialization
#include <iostream> // For error messages

NodeData::NodeData(int q_vel) : f(q_vel), f_eq(q_vel), f_new(q_vel) {}

void NodeData::compute_macroscopics(const Lattice& lattice) {
    if (!is_fluid) {
        rho = 0.0; // Or some other indicator for non-fluid nodes
        u = {0.0, 0.0};
        return;
    }

    rho = 0.0;
    u = {0.0, 0.0};
    for (int i = 0; i < lattice.get_Q(); ++i) {
        rho += f[i];
        u += lattice.get_c()[i] * f[i];
    }
    if (rho > REAL_EPSILON) { // Avoid division by zero or near-zero
        u /= rho;
    } else {
        // Handle potential instability or vacuum
        rho = REAL_EPSILON; // Reset to a small positive value? Or handle error.
        u = {0.0, 0.0};
        // Optionally set f to equilibrium for this small rho?
    }
}

void NodeData::initialize_equilibrium(const Lattice& lattice) {
     if (!EntropicEquilibrium::compute(*this, lattice)) {
         std::cerr << "Warning: Failed to compute initial equilibrium. Using polynomial approx." << std::endl;
         // Fallback to polynomial equilibrium (Eq. 5 in paper) if entropic fails
         Real u_sq = magnitude_sq(u);
         for (int i = 0; i < lattice.get_Q(); ++i) {
             Real cu = dot(lattice.get_c()[i], u);
             f_eq[i] = lattice.get_w()[i] * rho * (1.0 + cu / lattice.get_cs2() + 0.5 * (cu * cu) / (lattice.get_cs2() * lattice.get_cs2()) - 0.5 * u_sq / lattice.get_cs2());
         }
     }
     // Initialize f to f_eq
     f = f_eq;
}

Real NodeData::calculate_equilibrium(int k, const Lattice& lattice) const {
    // Calculate polynomial equilibrium for direction k
    Real u_sq = magnitude_sq(u);
    Real cu = dot(lattice.get_c()[k], u);
    Real cs2 = lattice.get_cs2();

    return lattice.get_w()[k] * rho * (1.0 + cu / cs2 + 0.5 * (cu * cu) / (cs2 * cs2) - 0.5 * u_sq / cs2);
}