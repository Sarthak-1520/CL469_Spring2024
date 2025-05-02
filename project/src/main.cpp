#include "simulator.h"
#include <iostream>
#include <cmath> // Required for std::sqrt

int main() {
    std::cout << "Starting Entropic Lattice Boltzmann Simulation (Flow Around Ellipse)..." << std::endl;

    // --- Simulation Parameters --- 
    int nx = 300;          // Grid resolution in X (make larger for wake)
    int ny = 100;          // Grid resolution in Y
    Real Re = 50.0;        // Reynolds number (e.g., 50-100 for interesting wake)

    // Ellipse parameters
    Real ellipse_cx = nx / 4.0;      // Center X (place upstream)
    Real ellipse_cy = ny / 2.0;      // Center Y
    Real ellipse_a = ny / 8.0;       // Semi-axis X
    Real ellipse_b = ny / 8.0;       // Semi-axis Y (make it a circle for simplicity first? Let's keep ellipse)

    // Inlet velocity based on low Mach number
    Real cs = 1.0 / std::sqrt(3.0);
    Real target_Ma = 0.05;
    Real U_in = target_Ma * cs; // Inlet velocity

    // Characteristic length L = 2 * ellipse_b (minor diameter)
    Real characteristic_length = 2.0 * ellipse_b;
    // viscosity = U * L / Re
    Real viscosity = U_in * characteristic_length / Re;

    int total_steps = 30000;   // Total simulation steps (might need more for stable wake)
    int output_freq = 500;    // Frequency of writing output files

    std::cout << "--- Ellipse Parameters ---" << std::endl;
    std::cout << "Center: (" << ellipse_cx << ", " << ellipse_cy << ")" << std::endl;
    std::cout << "Semi-axes: a=" << ellipse_a << ", b=" << ellipse_b << std::endl;
    std::cout << "Inlet U: " << U_in << " (Ma=" << U_in/cs << ")" << std::endl;
    std::cout << "Re = " << Re << ", Viscosity = " << viscosity << std::endl;
    std::cout << "-------------------------" << std::endl;

    // --- Create and Run Simulator ---
    try {
        // Pass ellipse parameters to the simulator constructor (needs modification)
        // For now, we'll handle geometry inside the simulator methods.
        Simulator sim(nx, ny, viscosity, U_in, total_steps, output_freq);
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