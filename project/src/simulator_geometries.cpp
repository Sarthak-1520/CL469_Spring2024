#include "simulator.h"
#include <iostream>
#include <cmath>

// ===== Elliptical Flow Setup =====
void Simulator::setup_ellipse_flow() {
    // Extract ellipse parameters from geometry parameters
    Real ellipse_cx = geom_param1 > 0.0 ? geom_param1 : nx / 4.0;
    Real ellipse_cy = geom_param2 > 0.0 ? geom_param2 : ny / 2.0;
    Real ellipse_a = geom_param3 > 0.0 ? geom_param3 : ny / 8.0;
    Real ellipse_b = geom_param4 > 0.0 ? geom_param4 : ny / 8.0;
    
    std::cout << "Setting up elliptical obstacle at (" << ellipse_cx << ", " << ellipse_cy << ") "
              << "with semi-axes a=" << ellipse_a << ", b=" << ellipse_b << std::endl;
    
    // Mark nodes inside the ellipse as non-fluid
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            Real term1 = ((static_cast<Real>(i) - ellipse_cx) / ellipse_a);
            Real term2 = ((static_cast<Real>(j) - ellipse_cy) / ellipse_b);
            if ((term1 * term1 + term2 * term2) <= 1.0) {
                grid[i][j].is_fluid = false;
            }
        }
    }

    // Mark top and bottom walls as non-fluid
    for (int i = 0; i < nx; ++i) {
        grid[i][0].is_fluid = false;      // Bottom wall
        grid[i][ny - 1].is_fluid = false; // Top wall
    }

    // Inlet (i=0) and Outlet (i=nx-1) remain fluid for now; handled by apply_all
    // Ensure corners are marked correctly based on wall logic
    grid[0][0].is_fluid = false;
    grid[nx-1][0].is_fluid = false;
    grid[0][ny-1].is_fluid = false;
    grid[nx-1][ny-1].is_fluid = false;
}

// ===== Lid-Driven Cavity Setup =====
void Simulator::setup_lid_driven_cavity() {
    std::cout << "Setting up lid-driven cavity..." << std::endl;
    
    // Mark all walls as non-fluid except the top (lid)
    for (int i = 0; i < nx; ++i) {
        grid[i][0].is_fluid = false;      // Bottom wall
        
        // Top wall (lid) remains fluid but will have special BC
        // We'll handle this in the boundary conditions
    }
    
    for (int j = 0; j < ny; ++j) {
        grid[0][j].is_fluid = false;      // Left wall
        grid[nx-1][j].is_fluid = false;   // Right wall
    }
    
    // Ensure corners are marked correctly
    grid[0][0].is_fluid = false;
    grid[nx-1][0].is_fluid = false;
    grid[0][ny-1].is_fluid = false;
    grid[nx-1][ny-1].is_fluid = false;
}

// ===== Channel with Obstacle Setup =====
void Simulator::setup_channel_obstacle() {
    // Extract obstacle parameters from geometry parameters
    Real obstacle_cx = geom_param1 > 0.0 ? geom_param1 : nx / 3.0;
    Real obstacle_cy = geom_param2 > 0.0 ? geom_param2 : ny / 2.0;
    Real obstacle_r = geom_param3 > 0.0 ? geom_param3 : ny / 8.0;
    
    std::cout << "Setting up circular obstacle at (" << obstacle_cx << ", " << obstacle_cy << ") "
              << "with radius r=" << obstacle_r << std::endl;
    
    // Mark nodes inside the circular obstacle as non-fluid
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            Real dx = static_cast<Real>(i) - obstacle_cx;
            Real dy = static_cast<Real>(j) - obstacle_cy;
            if ((dx * dx + dy * dy) <= obstacle_r * obstacle_r) {
                grid[i][j].is_fluid = false;
            }
        }
    }

    // Mark top and bottom walls as non-fluid
    for (int i = 0; i < nx; ++i) {
        grid[i][0].is_fluid = false;      // Bottom wall
        grid[i][ny - 1].is_fluid = false; // Top wall
    }

    // Inlet (i=0) and Outlet (i=nx-1) remain fluid for now; handled by apply_all
    // Ensure corners are marked correctly based on wall logic
    grid[0][0].is_fluid = false;
    grid[nx-1][0].is_fluid = false;
    grid[0][ny-1].is_fluid = false;
    grid[nx-1][ny-1].is_fluid = false;
}

// ===== Backward-Facing Step Setup =====
void Simulator::setup_backward_facing_step() {
    // Extract step parameters from geometry parameters
    int step_length = static_cast<int>(geom_param1 > 0.0 ? geom_param1 : nx / 6.0);
    int step_height = static_cast<int>(geom_param2 > 0.0 ? geom_param2 : ny / 2.0);
    
    std::cout << "Setting up backward-facing step with length=" << step_length 
              << " and height=" << step_height << std::endl;
    
    // Mark the step region as non-fluid
    for (int i = 0; i < step_length; ++i) {
        for (int j = 0; j < step_height; ++j) {
            grid[i][j].is_fluid = false;
        }
    }
    
    // Mark top and bottom walls as non-fluid
    for (int i = 0; i < nx; ++i) {
        grid[i][0].is_fluid = false;      // Bottom wall
        grid[i][ny - 1].is_fluid = false; // Top wall
    }
    
    // Inlet (i=0) and Outlet (i=nx-1) remain fluid for now; handled by apply_all
    // Ensure corners are marked correctly
    grid[0][0].is_fluid = false;
    grid[nx-1][0].is_fluid = false;
    grid[0][ny-1].is_fluid = false;
    grid[nx-1][ny-1].is_fluid = false;
}

// ===== Taylor-Green Vortex Setup =====
void Simulator::setup_taylor_green_vortex() {
    std::cout << "Setting up Taylor-Green vortex with periodic boundaries..." << std::endl;
    
    // All nodes are fluid for Taylor-Green vortex
    // Periodic boundary conditions will be applied
    
    // No additional setup needed as all nodes are already marked as fluid
}

// ===== Poiseuille Flow Setup =====
void Simulator::setup_poiseuille_flow() {
    std::cout << "Setting up Poiseuille flow channel..." << std::endl;
    
    // Mark top and bottom walls as non-fluid
    for (int i = 0; i < nx; ++i) {
        grid[i][0].is_fluid = false;      // Bottom wall
        grid[i][ny - 1].is_fluid = false; // Top wall
    }
    
    // Inlet (i=0) and Outlet (i=nx-1) remain fluid for now; handled by apply_all
    // Ensure corners are marked correctly
    grid[0][0].is_fluid = false;
    grid[nx-1][0].is_fluid = false;
    grid[0][ny-1].is_fluid = false;
    grid[nx-1][ny-1].is_fluid = false;
}
