#include "entropy.h"
#include <cmath>
#include <iostream> // For warnings

Real EntropyCalculator::calculate_H(const NodeData& node, const Lattice& lattice) {
    return calculate_H(node.f, lattice);
}

Real EntropyCalculator::calculate_H(const std::vector<Real>& f_vec, const Lattice& lattice) {
    Real H = 0.0;
    const auto& w = lattice.get_w();
    int Q = lattice.get_Q();

    for (int i = 0; i < Q; ++i) {
        if (f_vec[i] > REAL_EPSILON && w[i] > REAL_EPSILON) { // Check for positive f_i and w_i
            H += f_vec[i] * std::log(f_vec[i] / w[i]);
        } else if (f_vec[i] > 0.0 && f_vec[i] <= REAL_EPSILON) {
             // Treat very small positive f_i as contributing zero to avoid log(small number) issues
             // H += 0.0;
        } else if (f_vec[i] <= 0.0) {
            // This should ideally not happen in a stable simulation after collision
            // std::cerr << "Warning: Non-positive f_i encountered in H calculation: f[" << i << "] = " << f_vec[i] << std::endl;
             // Return a large value to indicate an invalid state?
             return REAL_MAX;
        }
        // If w[i] is zero or negative, the lattice definition is problematic
    }
    return H;
}