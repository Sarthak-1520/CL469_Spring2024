#include "entropic_equilibrium.h"
#include "numerical_solvers.h" // For solve_3x3
#include <cmath>
#include <iostream> // For debugging/warnings

// --- NewtonState Member Functions ---

void EntropicEquilibrium::NewtonState::update_residual_jacobian() {
    // Reset residual and jacobian
    residual = {0.0, 0.0, 0.0};
    jacobian = {}; // Zero initialize

    Real current_rho_eq = 0.0;
    Vector2D current_mom_eq = {0.0, 0.0};
    Real sum_w_c0_c0_expL = 0.0;
    Real sum_w_c0_c1_expL = 0.0;
    Real sum_w_c1_c1_expL = 0.0;

    const auto& c = lattice.get_c();
    const auto& w = lattice.get_w();
    int Q = lattice.get_Q();

    for (int i = 0; i < Q; ++i) {
        Real L = lambda[0] + lambda[1] * c[i].x + lambda[2] * c[i].y;
        Real expL = std::exp(L);
        Real w_expL = w[i] * expL;

        current_rho_eq += w_expL;
        current_mom_eq.x += w_expL * c[i].x;
        current_mom_eq.y += w_expL * c[i].y;

        // Jacobian terms (lower triangle, symmetric)
        sum_w_c0_c0_expL += w_expL * c[i].x * c[i].x;
        sum_w_c0_c1_expL += w_expL * c[i].x * c[i].y;
        sum_w_c1_c1_expL += w_expL * c[i].y * c[i].y;
    }

    // Residual vector R = [sum(f_eq)-rho, sum(cx*f_eq)-rho*ux, sum(cy*f_eq)-rho*uy]
    residual[0] = current_rho_eq - rho;
    residual[1] = current_mom_eq.x - rho * u.x;
    residual[2] = current_mom_eq.y - rho * u.y;

    // Jacobian matrix J = dR/dLambda
    jacobian[0][0] = current_rho_eq;
    jacobian[0][1] = current_mom_eq.x;
    jacobian[0][2] = current_mom_eq.y;
    jacobian[1][0] = jacobian[0][1]; // Symmetric
    jacobian[1][1] = sum_w_c0_c0_expL;
    jacobian[1][2] = sum_w_c0_c1_expL;
    jacobian[2][0] = jacobian[0][2]; // Symmetric
    jacobian[2][1] = jacobian[1][2]; // Symmetric
    jacobian[2][2] = sum_w_c1_c1_expL;
}

bool EntropicEquilibrium::NewtonState::solve_linear_system(std::array<Real, 3>& delta_lambda) {
    // Negate residual because we solve J * delta_lambda = -R
    std::array<Real, 3> neg_residual = {-residual[0], -residual[1], -residual[2]};
    return NumericalSolvers::solve_3x3(jacobian, neg_residual, delta_lambda);
}


// --- EntropicEquilibrium Static Member Function ---

bool EntropicEquilibrium::compute(NodeData& node, const Lattice& lattice, int max_iter, Real tol) {
    if (!node.is_fluid) return true; // Nothing to compute for non-fluid

    NewtonState state(node.rho, node.u, lattice);
    Real initial_residual_norm = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        state.update_residual_jacobian();

        Real current_residual_norm_sq = state.residual[0] * state.residual[0] +
                                        state.residual[1] * state.residual[1] +
                                        state.residual[2] * state.residual[2];
        Real current_residual_norm = std::sqrt(current_residual_norm_sq);

        if (iter == 0) {
            initial_residual_norm = current_residual_norm;
            if (initial_residual_norm < tol) break; // Already converged
        }

        // Check for convergence
        if (current_residual_norm < tol * (initial_residual_norm + tol)) { // Relative + absolute tolerance
             break; // Converged
        }

        // Solve for update step delta_lambda
        std::array<Real, 3> delta_lambda;
        if (!state.solve_linear_system(delta_lambda)) {
            std::cerr << "Warning: Jacobian is singular during f_eq calculation. Aborting update." << std::endl;
            return false; // Solver failed
        }

        // Update lambda (Lagrange multipliers)
        state.lambda[0] += delta_lambda[0];
        state.lambda[1] += delta_lambda[1];
        state.lambda[2] += delta_lambda[2];

        // Check for divergence (optional)
        if (iter == max_iter - 1) {
             std::cerr << "Warning: Newton solver for f_eq did not converge within " << max_iter << " iterations. Residual norm: " << current_residual_norm << std::endl;
             return false; // Failed to converge
        }
    }

    // Compute final f_eq using converged lambda values
    const auto& c = lattice.get_c();
    const auto& w = lattice.get_w();
    int Q = lattice.get_Q();
    for (int i = 0; i < Q; ++i) {
        Real L = state.lambda[0] + state.lambda[1] * c[i].x + state.lambda[2] * c[i].y;
        node.f_eq[i] = w[i] * std::exp(L);
        // Check for negative f_eq (shouldn't happen with exp)
        if (node.f_eq[i] < 0.0) {
             std::cerr << "Error: Negative f_eq[" << i << "] = " << node.f_eq[i] << " computed!" << std::endl;
             node.f_eq[i] = REAL_EPSILON; // Force positivity? Or return false.
             // return false;
        }
    }
    return true; // Success
}