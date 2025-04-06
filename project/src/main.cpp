#include "simulator.h"
#include <iostream>

int main() {
    std::cout << "Starting Entropic Lattice Boltzmann Simulation..." << std::endl;

    // --- Simulation Parameters (Lid-Driven Cavity Example) ---
    int nx = 100;          // Grid resolution in X
    int ny = 100;          // Grid resolution in Y
    Real Re = 1000000.0;       // Reynolds number

    // Characteristic length L = ny - 1 (height of cavity)
    // Characteristic velocity U = lid_velocity
    // Re = U * L / viscosity => viscosity = U * L / Re
    // We need to choose U (lid_velocity) such that Mach number Ma = U/cs is low (e.g., < 0.1)
    // cs = 1/sqrt(3) approx 0.577
    // Let Ma = 0.05 => U = Ma * cs = 0.05 / sqrt(3) approx 0.0288
    Real lid_velocity = 0.05 / std::sqrt(3.0); // Target Ma ~ 0.05
    Real characteristic_length = static_cast<Real>(ny - 1);
    Real viscosity = lid_velocity * characteristic_length / Re;

    int total_steps = 20000;   // Total simulation steps
    int output_freq = 500;    // Frequency of writing output files

    // --- Create and Run Simulator ---
    try {
        Simulator sim(nx, ny, viscosity, lid_velocity, total_steps, output_freq);
        sim.run();
    } catch (const std::exception& e) {
        std::cerr << "Error: Simulation failed with exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Error: Simulation failed with unknown exception." << std::endl;
        return 1;
    }


    std::cout << "Simulation completed successfully." << std::endl;
    return 0;
}