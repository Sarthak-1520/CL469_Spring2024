#pragma once
#include <iostream>
#include "node_data.h"
#include "common_types.h"
#include <vector>
#include <string>
#include <fstream>
#include <iomanip> // For setting precision

namespace VTKWriter {
    // Writes grid data (rho, u) to a VTK legacy file
    inline void write_vtk(const std::string& filename,
                          const std::vector<std::vector<NodeData>>& grid,
                          int nx, int ny)
    {
        std::ofstream vtk_file(filename);
        if (!vtk_file) {
            std::cerr << "Error: Could not open VTK file " << filename << std::endl;
            return;
        }

        vtk_file << "# vtk DataFile Version 3.0\n";
        vtk_file << "LBM Simulation Data\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET STRUCTURED_POINTS\n";
        vtk_file << "DIMENSIONS " << nx << " " << ny << " 1\n";
        vtk_file << "ORIGIN 0 0 0\n";
        vtk_file << "SPACING 1 1 1\n"; // Assuming dx=dy=dz=1
        vtk_file << "POINT_DATA " << nx * ny << "\n";

        // Write Density
        vtk_file << "SCALARS density double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                vtk_file << std::fixed << std::setprecision(10) << (grid[i][j].is_fluid ? grid[i][j].rho : 0.0) << "\n";
            }
        }

        // Write Velocity
        vtk_file << "VECTORS velocity double\n";
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                 vtk_file << std::fixed << std::setprecision(10)
                          << (grid[i][j].is_fluid ? grid[i][j].u.x : 0.0) << " "
                          << (grid[i][j].is_fluid ? grid[i][j].u.y : 0.0) << " 0.0\n"; // Z-component is 0 for 2D
            }
        }

        vtk_file.close();
    }

} // namespace VTKWriter