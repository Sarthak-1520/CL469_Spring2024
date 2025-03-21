#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip> 


#define PI 3.14159265358979323846

// D2Q9 model
const int Q = 9;  // Number of discrete velocities

// Discrete velocity vectors for D2Q9
const int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

// Opposite directions for bounce-back
const int opposite[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

// D2Q9 model weights
const double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
                    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

// Global parameters
int NX, NY;  // Lattice dimensions
int nsteps;  // Total number of time steps
int output_freq;  // Frequency of data output
double tau;  // Relaxation time
double omega;  // Relaxation parameter = 1/tau
double force_x, force_y;  // Applied force

class LatticeBoltzmann {
private:
    std::vector<std::vector<std::vector<double>>> f;  // Distribution function
    std::vector<std::vector<std::vector<double>>> f_new;  // New distribution function
    std::vector<std::vector<double>> rho;  // Density
    std::vector<std::vector<double>> ux;   // x-velocity
    std::vector<std::vector<double>> uy;   // y-velocity
    std::vector<std::vector<int>> obstacle;  // Obstacle marker (1 for wall, 0 for fluid)

public:
    LatticeBoltzmann() {
        // Initialize arrays
        f.resize(NX, std::vector<std::vector<double>>(NY, std::vector<double>(Q, 0.0)));
        f_new.resize(NX, std::vector<std::vector<double>>(NY, std::vector<double>(Q, 0.0)));
        rho.resize(NX, std::vector<double>(NY, 0.0));
        ux.resize(NX, std::vector<double>(NY, 0.0));
        uy.resize(NX, std::vector<double>(NY, 0.0));
        obstacle.resize(NX, std::vector<int>(NY, 0));
        
        // Initialize walls (top and bottom boundaries)
        for (int i = 0; i < NX; i++) {
            obstacle[i][0] = 1;  // Bottom wall
            obstacle[i][NY-1] = 1;  // Top wall
        }
        
        // Initialize distribution function
        initialize();
    }

    void initialize() {
        const double rho0 = 1.0;  // Initial density
        
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                // Skip walls
                if (obstacle[i][j] == 1) continue;
                
                // Initialize macroscopic values
                rho[i][j] = rho0;
                ux[i][j] = 0.0;
                uy[i][j] = 0.0;
                
                // Initialize distribution function with equilibrium
                for (int k = 0; k < Q; k++) {
                    double cu = cx[k] * ux[i][j] + cy[k] * uy[i][j];
                    double u2 = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
                    
                    f[i][j][k] = w[k] * rho[i][j] * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);
                    f_new[i][j][k] = f[i][j][k];
                }
            }
        }
    }

    void compute_macroscopic() {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                // Skip walls
                if (obstacle[i][j] == 1) continue;
                
                // Reset macroscopic values
                rho[i][j] = 0.0;
                ux[i][j] = 0.0;
                uy[i][j] = 0.0;
                
                // Calculate density and velocity
                for (int k = 0; k < Q; k++) {
                    rho[i][j] += f[i][j][k];
                    ux[i][j] += cx[k] * f[i][j][k];
                    uy[i][j] += cy[k] * f[i][j][k];
                }
                
                // Add half-step force effect (Guo's method)
                ux[i][j] = ux[i][j]/rho[i][j] + 0.5*force_x/rho[i][j];
                uy[i][j] = uy[i][j]/rho[i][j] + 0.5*force_y/rho[i][j];
            }
        }
    }

    void collision() {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                // Skip walls
                if (obstacle[i][j] == 1) continue;
                
                // Compute equilibrium distribution
                for (int k = 0; k < Q; k++) {
                    double cu = cx[k] * ux[i][j] + cy[k] * uy[i][j];
                    double u2 = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
                    double feq = w[k] * rho[i][j] * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);
                    
                    // Force term (Guo's method)
                    double force_term = 0.0;
                    if (k > 0) { // Skip rest particle
                        double cu = cx[k] * ux[i][j] + cy[k] * uy[i][j];
                        force_term = (1.0 - 0.5*omega) * w[k] * (
                            3.0 * ((cx[k] - ux[i][j])*force_x + (cy[k] - uy[i][j])*force_y) +
                            9.0 * (cx[k]*ux[i][j] + cy[k]*uy[i][j]) * (cx[k]*force_x + cy[k]*force_y)
                        );
                    }
                    
                    // BGK collision with force term
                    f_new[i][j][k] = f[i][j][k] - omega * (f[i][j][k] - feq) + force_term;
                }
            }
        }
    }

    void streaming() {
        // Create temporary copy for walls
        std::vector<std::vector<std::vector<double>>> f_temp = f_new;
        
        // Streaming step
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < Q; k++) {
                    int i_next = (i + cx[k] + NX) % NX;  // Periodic in x-direction
                    int j_next = j + cy[k];              // No periodicity in y-direction
                    
                    // Handle boundary conditions
                    if (j_next >= 0 && j_next < NY) {
                        if (obstacle[i_next][j_next] == 0) {
                            // Fluid node, regular streaming
                            f[i_next][j_next][k] = f_new[i][j][k];
                        } else {
                            // Wall node, bounce-back
                            f[i][j][opposite[k]] = f_new[i][j][k];
                        }
                    }
                }
            }
        }
    }

    void save_density(const std::string& filename) {
        std::ofstream file(filename);
        if (!file) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }
        
        // Use higher precision for output
        file << std::scientific << std::setprecision(10);
        
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                file << rho[i][j] << " ";
            }
            file << std::endl;
        }
        file.close();
        
        // Print debug info
        std::cout << "Saved density to " << filename << std::endl;
        
        // Find min and max density values
        double min_rho = 1e10, max_rho = -1e10;
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (rho[i][j] < min_rho) min_rho = rho[i][j];
                if (rho[i][j] > max_rho) max_rho = rho[i][j];
            }
        }
        std::cout << "  Density range: " << min_rho << " to " << max_rho << std::endl;
    }

    void save_velocity(const std::string& filename) {
        std::ofstream file(filename);
        if (!file) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }
        
        // Use higher precision for output
        file << std::scientific << std::setprecision(10);
        
        // Find min and max velocity values for debugging
        double min_ux = 1e10, max_ux = -1e10;
        double min_uy = 1e10, max_uy = -1e10;
        
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                // Save velocity components
                file << ux[i][j] << " " << uy[i][j] << " ";
                
                // Track min/max
                if (ux[i][j] < min_ux) min_ux = ux[i][j];
                if (ux[i][j] > max_ux) max_ux = ux[i][j];
                if (uy[i][j] < min_uy) min_uy = uy[i][j];
                if (uy[i][j] > max_uy) max_uy = uy[i][j];
            }
            file << std::endl;
        }
        file.close();
        
        // Print debug info
        std::cout << "Saved velocity to " << filename << std::endl;
        std::cout << "  Velocity range: ux = " << min_ux << " to " << max_ux 
                  << ", uy = " << min_uy << " to " << max_uy << std::endl;
        
        // Check if all velocities are 0 (which would indicate a problem)
        if (std::abs(min_ux) < 1e-10 && std::abs(max_ux) < 1e-10 && 
            std::abs(min_uy) < 1e-10 && std::abs(max_uy) < 1e-10) {
            std::cerr << "WARNING: All velocity values are close to zero!" << std::endl;
            std::cerr << "This may indicate a problem with the simulation." << std::endl;
        }
    }

    void save_velocity_profile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }
        
        // Use higher precision for output
        file << std::scientific << std::setprecision(10);
        
        // Take profile at the middle of the channel
        int i_mid = NX / 2;
        
        // Calculate analytical solution for comparison
        for (int j = 0; j < NY; j++) {
            double y = static_cast<double>(j);
            double h = static_cast<double>(NY - 1);
            double y_norm = y / h;  // Normalized coordinate
            
            // Analytical solution for Poiseuille flow
            double umax = force_x * h * h / (8 * (tau - 0.5) / 3.0);
            double u_analytical = umax * 4.0 * y_norm * (1.0 - y_norm);
            
            // For wall nodes, velocity should be 0
            if (obstacle[i_mid][j] == 1) {
                file << j << " " << 0.0 << " " << 0.0 << std::endl;
            } else {
                file << j << " " << ux[i_mid][j] << " " << u_analytical << std::endl;
            }
        }
        file.close();
        
        // Print debug info
        std::cout << "Saved velocity profile to " << filename << std::endl;
        std::cout << "  Profile taken at x = " << i_mid << std::endl;
        
        // Calculate and print max simulation and analytical velocities
        double max_sim_vel = 0.0;
        double max_ana_vel = 0.0;
        int j_max_sim = 0;
        
        for (int j = 0; j < NY; j++) {
            if (j == 0 || j == NY-1) continue; // Skip walls
            
            if (std::abs(ux[i_mid][j]) > max_sim_vel) {
                max_sim_vel = std::abs(ux[i_mid][j]);
                j_max_sim = j;
            }
            
            double y = static_cast<double>(j);
            double h = static_cast<double>(NY - 1);
            double y_norm = y / h;
            double umax = force_x * h * h / (8 * (tau - 0.5) / 3.0);
            double u_analytical = umax * 4.0 * y_norm * (1.0 - y_norm);
            
            if (std::abs(u_analytical) > max_ana_vel) {
                max_ana_vel = std::abs(u_analytical);
            }
        }
        
        std::cout << "  Max simulation velocity: " << max_sim_vel << " at j = " << j_max_sim << std::endl;
        std::cout << "  Max analytical velocity: " << max_ana_vel << std::endl;
        std::cout << "  Relative error: " << (max_sim_vel - max_ana_vel) / max_ana_vel * 100.0 << "%" << std::endl;
    }

    void run_simulation() {
        std::cout << "Starting simulation with parameters:" << std::endl;
        std::cout << "NX = " << NX << ", NY = " << NY << std::endl;
        std::cout << "tau = " << tau << ", omega = " << omega << std::endl;
        std::cout << "force_x = " << force_x << ", force_y = " << force_y << std::endl;
        std::cout << "nsteps = " << nsteps << ", output_freq = " << output_freq << std::endl;
        
        // Create output directory
        system("mkdir -p output");
        
        // Save initial state
        compute_macroscopic();
        std::string density_file = "output/rho_0.dat";
        std::string velocity_file = "output/vel_0.dat";
        std::string profile_file = "output/profile_0.dat";
        
        save_density(density_file);
        save_velocity(velocity_file);
        save_velocity_profile(profile_file);
        
        for (int step = 1; step <= nsteps; step++) {
            collision();
            streaming();
            compute_macroscopic();
            
            // Save results at specified intervals
            if (step % output_freq == 0) {
                std::cout << "Step " << step << "/" << nsteps << std::endl;
                
                std::string density_file = "output/rho_" + std::to_string(step) + ".dat";
                std::string velocity_file = "output/vel_" + std::to_string(step) + ".dat";
                std::string profile_file = "output/profile_" + std::to_string(step) + ".dat";
                
                save_density(density_file);
                save_velocity(velocity_file);
                save_velocity_profile(profile_file);
            }
        }
        
        // Save final results
        save_velocity_profile("output/final_profile.dat");
        std::cout << "Simulation completed successfully." << std::endl;
    }
};

void read_parameters() {
    std::ifstream file("parameters.dat");
    if (!file) {
        std::cerr << "Error: Could not open parameters.dat" << std::endl;
        exit(1);
    }
    
    file >> NX >> NY;
    file >> nsteps >> output_freq;
    file >> tau;
    file >> force_x >> force_y;
    
    omega = 1.0 / tau;
    
    std::cout << "Parameters loaded:" << std::endl;
    std::cout << "NX = " << NX << ", NY = " << NY << std::endl;
    std::cout << "nsteps = " << nsteps << ", output_freq = " << output_freq << std::endl;
    std::cout << "tau = " << tau << ", omega = " << omega << std::endl;
    std::cout << "force_x = " << force_x << ", force_y = " << force_y << std::endl;
    
    // Verify parameters for common issues
    if (NX <= 0 || NY <= 0) {
        std::cerr << "Error: Invalid grid dimensions!" << std::endl;
        exit(1);
    }
    
    if (tau <= 0.5) {
        std::cerr << "Warning: tau <= 0.5 can lead to numerical instability!" << std::endl;
    }
    
    if (force_x == 0 && force_y == 0) {
        std::cerr << "Warning: Both force components are zero. No flow will develop!" << std::endl;
    }
    
    file.close();
}

int main() {
    // Read parameters from file
    read_parameters();
    
    // Initialize and run simulation
    LatticeBoltzmann lbm;
    lbm.run_simulation();
    
    return 0;
} 