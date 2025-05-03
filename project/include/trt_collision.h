#pragma once
#include "node_data.h"
#include "lattice.h"

/**
 * @brief Two-Relaxation-Time (TRT) collision operator
 * 
 * The TRT collision operator separates the relaxation of symmetric (even) and antisymmetric (odd)
 * parts of the distribution functions, allowing for independent control of viscosity and numerical stability.
 * 
 * This implementation follows the approach described in:
 * Ginzburg, I., Verhaeghe, F., d'Humières, D. (2008). Two-relaxation-time Lattice Boltzmann scheme:
 * About parametrization, velocity, pressure and mixed boundary conditions.
 * Communications in Computational Physics, 3(2), 427-478.
 */
class TRTCollision {
public:
    /**
     * @brief Performs the TRT collision step for a single node
     * 
     * @param node The node to perform collision on
     * @param lattice The lattice definition
     * @param tau_plus Relaxation time for symmetric part (related to viscosity)
     * @param tau_minus Relaxation time for antisymmetric part (typically set using magic parameter)
     * @return void
     */
    static void collide(NodeData& node, const Lattice& lattice, Real tau_plus, Real tau_minus);

    /**
     * @brief Computes the magic parameter Lambda = (tau_plus - 0.5) * (tau_minus - 0.5)
     * 
     * @param tau_plus Relaxation time for symmetric part
     * @return Real The corresponding tau_minus value
     */
    static Real compute_tau_minus(Real tau_plus, Real magic_param = 1.0/4.0);
};
