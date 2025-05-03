#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <map>

// Structure to hold simulation parameters
struct SimParams {
    // Lattice Units
    double tau;
    double nu;
    int H; // Channel height in lattice units (Ny)
    int Lx; // Channel length in lattice units (Nx, will set to H)
    double g; // Body force (gravity) in x-direction
    double rho0; // Initial density
    int max_t_steps;
    int output_freq;
    std::string collision_operator; // "BGK" or "TRT"
    double tau_minus; // Free parameter for TRT (often 1.0 / tau_plus)

    // Physical Parameters (for reference)
    double Re_p;
    double nu_p;
    double H_p;
    double g_p;
    double rho_p;

    // Calculated Lattice parameters
    double cs2; // Speed of sound squared
    double omega; // Relaxation frequency (1/tau) for BGK
    double omega_plus; // Symmetric relaxation freq for TRT (1/tau_plus)
    double omega_minus; // Anti-symmetric relaxation freq for TRT (1/tau_minus)
    double force_x; // Lattice force component fx = g * rho0 (assuming constant density for force calc)
    double force_y; // Lattice force component fy = 0

    // Derived physical scales (for non-dimensionalization)
    double um_analytical; // Max analytical velocity in lattice units
    double F_analytical_nondim; // Analytical non-dim force = 0.5
};

// Function to read parameters from a file
SimParams readParameters(const std::string& filename);

// Helper to parse the parameter file
std::map<std::string, std::string> parseParamFile(const std::string& filename);

#endif // PARAMETERS_H 