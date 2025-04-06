#pragma once
#include "node_data.h"
#include "lattice.h"

class EntropicEquilibrium {
public:
    // Calculates f_eq based on rho, u using entropy minimization (Newton-Raphson)
    // Stores the result in node.f_eq
    // Returns true on success, false on failure (e.g., solver divergence)
    static bool compute(NodeData& node, const Lattice& lattice, int max_iter = 100, Real tol = 1e-7);

private:
    // Structure to hold state for Newton-Raphson solver
    struct NewtonState {
        Vector2D u;
        Real rho;
        const Lattice& lattice;
        std::array<Real, 3> lambda = {0.0, 0.0, 0.0}; // lambda_0, lambda_x, lambda_y
        std::array<Real, 3> residual = {0.0, 0.0, 0.0};
        std::array<std::array<Real, 3>, 3> jacobian = {}; // 3x3 matrix

        NewtonState(Real r, const Vector2D& vel, const Lattice& lat) : u(vel), rho(r), lattice(lat) {}

        void update_residual_jacobian();
        bool solve_linear_system(std::array<Real, 3>& delta_lambda);
    };
};