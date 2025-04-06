#pragma once
#include "node_data.h"
#include "common_types.h"
#include <vector>
#include <string>
#include <fstream>
#include <iomanip> // For setting precision
#include <iostream> // For std::cerr

namespace DataWriter {
    // Writes grid data (rho, u) to a simple text .dat file
    inline void write_dat(const std::string& filename,
                          const std::vector<std::vector<NodeData>>& grid,
                          int nx, int ny)
    {
        std::ofstream dat_file(filename);
        if (!dat_file) {
            std::cerr << "Error: Could not open DAT file " << filename << std::endl;
            return;
        }

        // Write header: NX NY
        dat_file << "# NX NY\n";
        dat_file << nx << " " << ny << "\n";

        // Write data header
        dat_file << "# i j rho u.x u.y is_fluid\n";

        // Write data for each point
        dat_file << std::fixed << std::setprecision(10);
        for (int j = 0; j < ny; ++j) { // Iterate y first for row-major like access if reading into numpy
            for (int i = 0; i < nx; ++i) {
                 dat_file << i << " " << j << " "
                          << grid[i][j].rho << " "
                          << grid[i][j].u.x << " "
                          << grid[i][j].u.y << " "
                          << (grid[i][j].is_fluid ? 1 : 0) // Add is_fluid flag
                          << "\n";
            }
        }

        dat_file.close();
    }

} // namespace DataWriter 