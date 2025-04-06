#include "lattice.h"

Lattice::Lattice() : cs2(1.0 / 3.0) {
    // Velocities for D2Q9
    c.resize(Q);
    c[0] = {0, 0};
    c[1] = {1, 0};  c[2] = {0, 1};  c[3] = {-1, 0}; c[4] = {0, -1};
    c[5] = {1, 1};  c[6] = {-1, 1}; c[7] = {-1, -1};c[8] = {1, -1};

    // Weights for D2Q9
    w.resize(Q);
    w[0] = 4.0 / 9.0;
    w[1] = w[2] = w[3] = w[4] = 1.0 / 9.0;
    w[5] = w[6] = w[7] = w[8] = 1.0 / 36.0;

    // Opposite directions
    opp[0] = 0; opp[1] = 3; opp[2] = 4; opp[3] = 1; opp[4] = 2;
    opp[5] = 7; opp[6] = 8; opp[7] = 5; opp[8] = 6;
}

int Lattice::opposite(int i) const {
    // Basic bounds check
    if (i >= 0 && i < Q) {
        return opp[static_cast<size_t>(i)];
    }
    // Handle error case, though ideally this shouldn't be reached with valid input
    return -1; // Or throw an exception
}