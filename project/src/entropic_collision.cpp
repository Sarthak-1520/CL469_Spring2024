#include "entropic_collision.h"
#include "numerical_solvers.h" // For Brent's method
#include <vector>
#include <cmath>
#include <iostream> // For warnings

// --- AlphaObjective Member Functions ---

EntropicCollision::AlphaObjective::AlphaObjective(
    const std::vector<Real>& f_pre_in, const std::vector<Real>& f_eq_in,
    const Lattice& lattice_in, Real beta_in)
    : f_pre(f_pre_in), f_eq(f_eq_in), lattice(lattice_in), beta(beta_in), Q(lattice_in.get_Q())
{
    H_pre = EntropyCalculator::calculate_H(f_pre, lattice);
    // Handle case where initial state is already invalid
    if (H_pre >= REAL_MAX) {
         std::cerr << "Warning: Initial H calculation failed in AlphaObjective." << std::endl;
    }
}

Real EntropicCollision::AlphaObjective::operator()(Real alpha) const {
     if (H_pre >= REAL_MAX) return REAL_MAX; // Initial state was bad

    std::vector<Real> f_post(Q);
    bool non_positive_f = false;
    for (int i = 0; i < Q; ++i) {
        // f_post = f_pre + alpha * beta * (f_eq - f_pre)
        // Let's redefine the collision as: f_post = f_pre + relaxation_param * (f_eq - f_pre)
        // where relaxation_param = alpha * beta. We solve for alpha using H constraint and apply it.
        f_post[i] = f_pre[i] + alpha * beta * (f_eq[i] - f_pre[i]); // Use alpha*beta for trial state

        if (f_post[i] <= 0.0) {
            non_positive_f = true;
            break; // No need to continue if one f_i is invalid
        }
    }

    if (non_positive_f) {
        // Return a large positive value to guide the solver away from non-physical alpha
        return REAL_MAX;
    }

    Real H_post = EntropyCalculator::calculate_H(f_post, lattice);
    if (H_post >= REAL_MAX) {
         // Calculation failed for this alpha
         return REAL_MAX;
    }

    // We want H(f_post) - H(f_pre) = 0
    return H_post - H_pre;
}


// --- EntropicCollision Static Member Function ---

Real EntropicCollision::collide(NodeData& node, const Lattice& lattice, Real beta, Real alpha_tol, int alpha_max_iter) {
    if (!node.is_fluid) return 1.0; // No collision for non-fluid

    // --- Pre-Collision Positivity Clamp (Defensive) ---
    // Ensure f is positive before calculating H_pre, as non-positivity might come from streaming/BCs
    for (int i = 0; i < lattice.get_Q(); ++i) {
        if (node.f[i] <= 0.0) {
            // Uncomment the warning if you want to know when this pre-clamp is triggered
            // std::cerr << "Warning: Clamping non-positive f[" << i << "] = " << node.f[i] << " BEFORE H_pre calculation." << std::endl;
            node.f[i] = REAL_EPSILON;
        }
    }
    // --- End Pre-Collision Clamp ---

    // 1. Compute Entropic Equilibrium f_eq (should be done before calling collide)
    //    Assuming node.f_eq is already computed and valid.

    // 2. Define the objective function G(alpha) = H(f_post) - H(f_pre)
    AlphaObjective objective(node.f, node.f_eq, lattice, beta); // Now node.f is guaranteed positive

    // Check if initial H calculation failed (Could still fail if f_eq is bad, but less likely)
    if (objective.H_pre >= REAL_MAX) {
         std::cerr << "Error: Cannot perform entropic collision - H_pre calculation failed even after clamping f_pre." << std::endl;
         return 1.0;
    }

    // 3. Find the root alpha >= 1 using Brent's method
    //    Search range [1.0, alpha_max]. alpha_max=2 corresponds to standard LBGK limit.
    //    We need a slightly larger range to allow for the entropic adjustment.
    Real alpha_min_search = 1.0;
    Real alpha_max_search = 2.0; // Reduce search range slightly from 2.1 to 2.0

    // Evaluate function at bounds to ensure a root might exist (sign change)
    Real G_at_min = objective(alpha_min_search);
    Real G_at_max = objective(alpha_max_search);

    // Check if bounds are valid
     if (G_at_min >= REAL_MAX || G_at_max >= REAL_MAX) {
         std::cerr << "Warning: Objective function invalid at search bounds for alpha. Defaulting alpha=1.0" << std::endl;
         // Don't modify f, just return 1.0
         return 1.0;
     }

    std::optional<Real> alpha_opt;
    // If G(1) is already very close to zero, alpha=1 is the solution
    if (std::abs(G_at_min) < alpha_tol) {
        alpha_opt = 1.0;
    }
    // If G(1) > 0 (entropy increases at alpha=1, which shouldn't happen ideally)
    // or if G(1) and G(max) have the same sign, the root finding might fail or alpha=1 is the minimum entropy point.
    else if (G_at_min * G_at_max >= 0) {
         // No guaranteed root in [1, max]. This might happen if f is already very close to f_eq.
         // Or if the minimum entropy is exactly at alpha=1.
         // std::cerr << "Warning: No sign change for alpha objective in [" << alpha_min_search << ", " << alpha_max_search << "]. G(1)=" << G_at_min << ", G(max)=" << G_at_max << ". Defaulting alpha=1.0" << std::endl;
         alpha_opt = 1.0; // Default to alpha=1 if no root found or needed
    } else {
        // Proceed with Brent's method
        alpha_opt = NumericalSolvers::brent(objective, alpha_min_search, alpha_max_search, alpha_tol, alpha_max_iter);
    }


    Real alpha = 1.0; // Default value if solver fails
    if (alpha_opt) {
        alpha = *alpha_opt;
        // Ensure alpha is at least 1.0 due to potential numerical inaccuracies
        if (alpha < 1.0) {
            // std::cerr << "Warning: Brent solver returned alpha < 1.0 (" << alpha << "). Clamping to 1.0." << std::endl;
            alpha = 1.0;
        }
    } else {
        std::cerr << "Warning: Brent solver failed to find alpha. Defaulting alpha=1.0" << std::endl;
        alpha = 1.0;
    }

    // 4. Apply the collision step using the found alpha
    //    f_post = f_pre + alpha * beta * (f_eq - f_pre)  <- Original interpretation
    //    f_post = f_pre + alpha * (f_eq - f_pre)         <- Interpretation based on Eq 54 structure
    // Let's use the relaxation parameter alpha_eff = alpha * beta
    Real alpha_eff = alpha * beta;
    for (int i = 0; i < lattice.get_Q(); ++i) {
        node.f[i] = node.f[i] + alpha_eff * (node.f_eq[i] - node.f[i]);
        // Ensure positivity after collision (redundant if pre-clamp exists AND alpha solver works, but keep for safety)
        if (node.f[i] <= 0.0) { // Use <= 0.0 for safety
             // Uncomment the warning if desired for debugging, but clamp regardless
             // std::cerr << "Warning: Non-positive f[" << i << "] = " << node.f[i] << " after collision with alpha=" << alpha << ", beta=" << beta << ". Clamping." << std::endl;
             node.f[i] = REAL_EPSILON; // Clamp to small positive value
        }
    }

    return alpha; // Return the calculated alpha for potential analysis
}