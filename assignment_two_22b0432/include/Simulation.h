#ifndef SIMULATION_H
#define SIMULATION_H

#include "Parameters.h"
#include "Lattice.h"
#include <vector>
#include <string>
#include <fstream>

class Simulation {
public:
    Simulation(const SimParams& p);
    ~Simulation() = default;

    void run();

private:
    // Initialization
    void initialize();

    // Main steps
    void updateMacroscopic();
    void collide();
    void streamAndBC();

    // Collision operators
    void collideBGK();
    void collideTRT();

    // Equilibrium distribution
    void equilibrium(double rho_node, double ux_node, double uy_node, double* feq_node);

    // Force term calculation (Guo et al.)
    void calculateForceTermBGK(int i, double rho_node, double ux_node, double uy_node, double& Fi);
    void calculateForceTermTRT(int i, double ux_node, double uy_node, double& Fi_plus, double& Fi_minus);

    // Boundary conditions
    void applyBounceBack(int x, int y); // Applied during streaming

    // Analysis & Output
    void calculateForces();
    void calculateViscousDissipation();
    void writeData(int timestep);
    void writeResultsSummary();

    // Force calculation methods
    void calculateForcesMomentumExchange(double& F_bottom_x, double& F_bottom_y, double& F_top_x, double& F_top_y);
    void calculateStressTensor(std::vector<double>& tau_xx, std::vector<double>& tau_xy, std::vector<double>& tau_yy);
    void calculateForcesStressIntegration(double& F_bottom_x, double& F_bottom_y, double& F_top_x, double& F_top_y);
    void calculateForcesFiniteDifference(double& F_bottom_x, double& F_top_x);

    // Helper functions
    inline int idx(int x, int y, int i) const { return (y * Nx + x) * LBM::Q + i; }
    inline int idx_macro(int x, int y) const { return y * Nx + x; }
    double getVelocityChangeNorm(); // For checking convergence

    // Member variables
    SimParams params; // Simulation parameters
    int Nx, Ny; // Grid dimensions

    std::vector<double> f;      // Current populations f_i(x, y, t)
    std::vector<double> f_new;  // New populations f_i(x, y, t+dt)
    std::vector<double> rho;    // Density rho(x, y)
    std::vector<double> ux;     // Velocity u_x(x, y)
    std::vector<double> uy;     // Velocity u_y(x, y)
    std::vector<double> ux_old; // Velocity u_x(x, y) at previous step for convergence check

    // Results storage (accumulated or final)
    double F_bottom_x_me = 0.0, F_bottom_y_me = 0.0;
    double F_top_x_me = 0.0, F_top_y_me = 0.0;
    double F_bottom_x_si = 0.0, F_bottom_y_si = 0.0;
    double F_top_x_si = 0.0, F_top_y_si = 0.0;
    double F_bottom_x_fd = 0.0, F_top_x_fd = 0.0;
    double total_viscous_dissipation = 0.0;

    // File streams for output
    std::ofstream velocity_file;
    std::ofstream force_file;
    std::ofstream convergence_file;

};

#endif // SIMULATION_H 