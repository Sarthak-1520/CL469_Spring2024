#ifndef LATTICE_H
#define LATTICE_H

#include <vector>

namespace LBM {
    // D2Q9 Lattice constants
    const int Q = 9;
    const int D = 2;

    // Lattice velocities (c_i)
    const int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // Lattice weights (w_i)
    const double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
                         1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // Opposite direction index (opp[i] = index of -c_i)
    const int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Speed of sound squared
    const double cs2 = 1.0 / 3.0;

} // namespace LBM

#endif // LATTICE_H 