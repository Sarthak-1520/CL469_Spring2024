#pragma once
#include "node_data.h"
#include "lattice.h"
#include "entropy.h"

class EntropicCollision {
public:
    // Performs the entropic collision step for a single node
    // Updates node.f based on node.f (pre-collision) and node.f_eq
    // Returns the calculated alpha value (or 1.0 if solver fails)
    static Real collide(NodeData& node, const Lattice& lattice, Real beta /* dt / (2*tau + dt) */, Real alpha_tol = 1e-9, int alpha_max_iter = 50);

private:
    // Function defining the root-finding target for alpha (Eq. 54)
    // G(alpha) = H(f_pre + alpha * beta * (f_eq - f_pre)) - H(f_pre) = 0
    struct AlphaObjective {
        const std::vector<Real>& f_pre;
        const std::vector<Real>& f_eq;
        const Lattice& lattice;
        Real beta;
        Real H_pre;
        int Q;

        AlphaObjective(const std::vector<Real>& f_pre_in, const std::vector<Real>& f_eq_in,
                       const Lattice& lattice_in, Real beta_in);

        // Evaluate the objective function G(alpha)
        // Returns G(alpha) or REAL_MAX if f becomes non-positive
        Real operator()(Real alpha) const;
    };
};