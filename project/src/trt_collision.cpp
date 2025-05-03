#include "trt_collision.h"
#include <vector>
#include <cmath>
#include <iostream>

void TRTCollision::collide(NodeData& node, const Lattice& lattice, Real tau_plus, Real tau_minus) {
    if (!node.is_fluid) return; // No collision for non-fluid nodes

    const int Q = lattice.get_Q();

    // Ensure f_eq is computed before calling this function
    // We assume node.f_eq contains the equilibrium distribution

    // Temporary storage for post-collision distributions
    std::vector<Real> f_post(Q);

    // Compute symmetric and antisymmetric parts and apply relaxation
    for (int i = 0; i < Q; ++i) {
        int i_opp = lattice.opposite(i);

        // Symmetric part (even): f_i^+ = (f_i + f_{opp(i)}) / 2
        Real f_plus = 0.5 * (node.f[i] + node.f[i_opp]);
        Real f_eq_plus = 0.5 * (node.f_eq[i] + node.f_eq[i_opp]);

        // Antisymmetric part (odd): f_i^- = (f_i - f_{opp(i)}) / 2
        Real f_minus = 0.5 * (node.f[i] - node.f[i_opp]);
        Real f_eq_minus = 0.5 * (node.f_eq[i] - node.f_eq[i_opp]);

        // Apply TRT collision: relax symmetric and antisymmetric parts separately
        Real f_plus_post = f_plus - (f_plus - f_eq_plus) / tau_plus;
        Real f_minus_post = f_minus - (f_minus - f_eq_minus) / tau_minus;

        // Reconstruct post-collision distribution: f_i = f_i^+ + f_i^-
        f_post[i] = f_plus_post + f_minus_post;

        // Ensure positivity (defensive programming)
        if (f_post[i] <= 0.0) {
            // std::cerr << "Warning: Non-positive f_post[" << i << "] = " << f_post[i] << " in TRT collision. Clamping." << std::endl;
            f_post[i] = REAL_EPSILON;
        }
    }

    // Update node with post-collision distributions
    node.f = f_post;
}

Real TRTCollision::compute_tau_minus(Real tau_plus, Real magic_param) {
    // The "magic parameter" Lambda = (tau_plus - 0.5) * (tau_minus - 0.5)
    // Solving for tau_minus: tau_minus = magic_param / (tau_plus - 0.5) + 0.5
    return magic_param / (tau_plus - 0.5) + 0.5;
}
