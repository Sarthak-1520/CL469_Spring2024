#pragma once
#include "node_data.h"
#include "lattice.h"

class EntropyCalculator {
public:
    // Calculates H = sum(f_i * log(f_i / w_i)) for a node
    static Real calculate_H(const NodeData& node, const Lattice& lattice);
    // Calculates H for a given distribution vector f
    static Real calculate_H(const std::vector<Real>& f_vec, const Lattice& lattice);
};