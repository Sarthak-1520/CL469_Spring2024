#pragma once
#include "common_types.h"
#include <vector>
#include <array>

class Lattice {
public:
    // Constructor for D2Q9
    Lattice();

    int get_D() const { return D; }
    int get_Q() const { return Q; }
    const std::vector<Vector2D>& get_c() const { return c; } // Discrete velocities
    const std::vector<Real>& get_w() const { return w; }     // Weights
    Real get_cs2() const { return cs2; }                     // Speed of sound squared
    int opposite(int i) const;                               // Opposite direction index

private:
    const int D = 2; // Dimension
    const int Q = 9; // Number of velocities
    std::vector<Vector2D> c; // Discrete velocity vectors (e_i or c_i in paper)
    std::vector<Real> w;     // Lattice weights
    std::array<int, 9> opp;  // Opposite direction indices
    Real cs2;                // Speed of sound squared (1/3)
};