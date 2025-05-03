#include "Simulation.h"
#include "Parameters.h"
#include <iostream>
#include <string>
#include <algorithm> // Required for std::transform
#include <cctype>    // Required for ::tolower

// Helper function to convert string to lower case
std::string toLowerMain(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

int main() {
    std::cout << "Starting LBM Poiseuille Flow Simulation..." << std::endl;

    try {
        // Read parameters
        std::string param_file = "parameters.dat";
        SimParams params = readParameters(param_file);
        std::cout << "Parameters loaded from " << param_file << std::endl;
        std::cout << "  Grid: " << params.Lx << " x " << params.H << std::endl;
        std::cout << "  tau = " << params.tau << ", nu = " << params.nu << std::endl;
        std::cout << "  g = " << params.g << ", rho0 = " << params.rho0 << std::endl;
        std::cout << "  Max steps = " << params.max_t_steps << ", Output freq = " << params.output_freq << std::endl;
        std::cout << "  Collision Operator: " << params.collision_operator << std::endl;
        if (toLowerMain(params.collision_operator) == "trt") {
            std::cout << "    tau_minus = " << params.tau_minus << std::endl;
        }

        // Create and run simulation
        Simulation sim(params);
        sim.run();

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Simulation finished successfully." << std::endl;
    return 0;
} 