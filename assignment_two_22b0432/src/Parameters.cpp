#include "Parameters.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm> // for std::transform
#include <cctype>    // for ::tolower

// Helper function to trim whitespace from a string
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (std::string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

// Helper function to convert string to lower case
std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

// Parses the parameter file into a map
std::map<std::string, std::string> parseParamFile(const std::string& filename) {
    std::map<std::string, std::string> paramsMap;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Error opening parameter file: " + filename);
    }

    while (std::getline(file, line)) {
        // Remove comments (lines starting with # or text after #)
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        // Trim whitespace
        line = trim(line);

        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Find the equals sign
        size_t equalsPos = line.find('=');
        if (equalsPos == std::string::npos) {
            std::cerr << "Warning: Skipping invalid line in parameter file: " << line << std::endl;
            continue;
        }

        std::string key = trim(line.substr(0, equalsPos));
        std::string value = trim(line.substr(equalsPos + 1));

        if (key.empty() || value.empty()) {
            std::cerr << "Warning: Skipping invalid key/value pair in parameter file: " << line << std::endl;
            continue;
        }

        paramsMap[key] = value;
    }

    file.close();
    return paramsMap;
}

// Reads parameters and populates the SimParams struct
SimParams readParameters(const std::string& filename) {
    SimParams p;
    std::map<std::string, std::string> paramsMap = parseParamFile(filename);

    auto getDouble = [&](const std::string& key) {
        if (paramsMap.count(key)) return std::stod(paramsMap.at(key));
        throw std::runtime_error("Missing parameter in file: " + key);
    };

    auto getInt = [&](const std::string& key) {
        if (paramsMap.count(key)) return std::stoi(paramsMap.at(key));
        throw std::runtime_error("Missing parameter in file: " + key);
    };

    auto getString = [&](const std::string& key, const std::string& defaultValue = "") {
        if (paramsMap.count(key)) return paramsMap.at(key);
        if (!defaultValue.empty()) return defaultValue;
        throw std::runtime_error("Missing parameter in file: " + key);
    };

    try {
        // Read Lattice Unit parameters
        p.tau = getDouble("tau");
        p.nu = getDouble("nu");
        p.H = getInt("H");
        p.Lx = p.H; // Set Lx = H as default based on assignment figure
        p.g = getDouble("g");
        p.rho0 = getDouble("rho0");
        p.max_t_steps = getInt("max_t_steps");
        p.output_freq = getInt("output_freq");

        // Optional parameters with defaults
        p.collision_operator = getString("collision_operator", "BGK");
        p.tau_minus = getDouble("tau_minus"); // Default needed? Set later based on tau_plus?
                                              // For now, assume it's provided if TRT is used.

        // Read Physical Parameters
        p.Re_p = getDouble("Re_p");
        p.nu_p = getDouble("nu_p");
        p.H_p = getDouble("H_p");
        p.g_p = getDouble("g_p");
        p.rho_p = getDouble("rho_p");

        // Calculate derived parameters
        p.cs2 = 1.0 / 3.0;
        p.omega = 1.0 / p.tau;
        // TRT omegas require tau_plus, which depends on nu. Set tau_plus = tau for now.
        // Proper TRT setup might require adjusting input params or how they're derived.
        double tau_plus = p.nu / p.cs2 + 0.5; // As per assignment text for consistency if needed
        p.omega_plus = 1.0 / tau_plus;
        if (toLower(p.collision_operator) == "trt") {
             if (!paramsMap.count("tau_minus")) {
                 // Default choice for tau_minus if not specified (Ginzburg et al. magic parameter)
                 // lambda = (1/tau_minus - 0.5)*(1/tau_plus - 0.5) = 1/4
                 // (2/omega_minus - 0.5)*(2/omega_plus - 0.5) = 1/4
                 // Let's set omega_minus = 1.0 for simplicity if not given, which implies tau_minus=1.0
                 std::cout << "Warning: tau_minus not specified for TRT, using default tau_minus = 1.0 (omega_minus = 1.0)" << std::endl;
                 p.tau_minus = 1.0;
             }
             p.omega_minus = 1.0 / p.tau_minus;
         } else {
            p.tau_minus = 0.0; // Not used for BGK
            p.omega_minus = 0.0; // Not used for BGK
         }


        p.force_x = p.g * p.rho0; // Force per unit volume F = rho * g
        p.force_y = 0.0;

        // Calculate analytical max velocity (Eq. 12 relation)
        // um = g * H^2 / (8 * nu)
        // Here H is the *width* (number of intervals), so H_lattice = Ny-1 = p.H - 1
        double H_lattice_width = static_cast<double>(p.H - 1);
        p.um_analytical = p.g * H_lattice_width * H_lattice_width / (8.0 * p.nu);

        // Analytical non-dimensional force
        p.F_analytical_nondim = 0.5;

        // Validate parameters
        if (p.tau <= 0.5 && toLower(p.collision_operator) == "bgk") {
             std::cerr << "Warning: tau <= 0.5, potential instability for BGK." << std::endl;
        }
         if (p.H <= 2) {
            throw std::runtime_error("Channel height H must be > 2.");
        }

    } catch (const std::exception& e) {
        std::cerr << "Error parsing parameters: " << e.what() << std::endl;
        throw;
    }

    return p;
} 