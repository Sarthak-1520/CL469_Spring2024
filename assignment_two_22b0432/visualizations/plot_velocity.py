import numpy as np
import matplotlib.pyplot as plt
import os

# --- Simulation Parameters (should match parameters.dat used for the run) ---
# Need H, nu, g to calculate analytical profile and scales
# TODO: Read these from parameters.dat or results_summary.txt for robustness
H = 11      # Channel height in lattice units (Ny)
nu = 0.1    # Kinematic viscosity (lattice units)
g = 9.203e-6 # Gravity (lattice units)
rho0 = 1.0  # Reference density
# --------------------------------------------------------------------------

DATA_DIR = 'data/' # Use path relative to assignment_two root
FIG_DIR = 'figures/'
VELOCITY_FILE = os.path.join(DATA_DIR, 'velocity_profile.dat')
OUTPUT_FIG = os.path.join(FIG_DIR, 'velocity_profile.png')

# Ensure figure directory exists
os.makedirs(FIG_DIR, exist_ok=True)

# Analytical solution (Eq. 5 adapted to lattice units)
# ux(y) = (g / (2 * nu)) * y * (H_width - y)
# y ranges from 0 to H_width (inclusive lattice nodes)
# H_width = Ny - 1 = H - 1 (number of intervals)
H_width = H - 1

# Calculate analytical max velocity (at y = H_width / 2)
umax_analytical = (g / (8.0 * nu)) * H_width**2
print(f"Analytical max velocity (lattice units): {umax_analytical:.4e}")

def analytical_velocity(y_lattice, H_width_lattice, g_lattice, nu_lattice):
    """Calculates analytical Poiseuille velocity profile in lattice units."""
    # Ensure y is within bounds [0, H_width]
    y = np.clip(y_lattice, 0, H_width_lattice)
    return (g_lattice / (2.0 * nu_lattice)) * y * (H_width_lattice - y)

# Load simulation data
try:
    data = np.loadtxt(VELOCITY_FILE)
except OSError:
    print(f"Error: Could not find or read {VELOCITY_FILE}")
    print("Please run the simulation first.")
    exit()
except ValueError:
    print(f"Error: Could not parse data in {VELOCITY_FILE}. Is it complete?")
    exit()

if data.size == 0:
    print(f"Error: {VELOCITY_FILE} is empty.")
    exit()

# Extract data for the last timestep
last_timestep = data[-1, 0]
print(f"Plotting data for timestep: {int(last_timestep)}")
final_data = data[data[:, 0] == last_timestep]

y_sim = final_data[:, 1]
ux_sim = final_data[:, 2]
uy_sim = final_data[:, 3]

# Check simulation results consistency
max_ux_sim = np.max(ux_sim)
max_uy_sim_abs = np.max(np.abs(uy_sim))
print(f"Max simulated Ux: {max_ux_sim:.4e}")
print(f"Max absolute simulated Uy: {max_uy_sim_abs:.4e} (should be close to 0)")

# Create analytical profile for comparison
y_analytical_lattice = np.linspace(0, H_width, 100)
ux_analytical = analytical_velocity(y_analytical_lattice, H_width, g, nu)

# Non-dimensionalize
y_sim_nd = y_sim / H_width
ux_sim_nd = ux_sim / umax_analytical

y_analytical_nd = y_analytical_lattice / H_width
ux_analytical_nd = ux_analytical / umax_analytical

# Plotting Velocity
plt.figure(figsize=(8, 6))
plt.plot(ux_analytical_nd, y_analytical_nd, 'r-', label='Analytical Solution (Eq. 5)')
plt.plot(ux_sim_nd, y_sim_nd, 'bo', markerfacecolor='none', markersize=6, label=f'LBM Simulation (t={int(last_timestep)})')

plt.xlabel('$u_x / u_{max,analytical}$')
plt.ylabel('$y / H$')
plt.title('Steady-State Velocity Profile Comparison (Ux)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.xlim(left=-0.1) # Ensure origin is visible

plt.savefig(OUTPUT_FIG)
print(f"Velocity profile plot saved to {OUTPUT_FIG}")
plt.close() # Optional: close plot window if showing automatically

# Plotting Uy
OUTPUT_FIG_UY = os.path.join(FIG_DIR, 'uy_profile.png')
plt.figure(figsize=(8, 6))
plt.plot(uy_sim, y_sim, 'go-', markerfacecolor='none', markersize=6, label=f'LBM Simulation (t={int(last_timestep)})')
plt.xlabel('$u_y$ (Lattice Units)')
plt.ylabel('$y$ (Lattice Units)')
plt.title('Steady-State Velocity Profile (Uy)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

plt.savefig(OUTPUT_FIG_UY)
print(f"Uy profile plot saved to {OUTPUT_FIG_UY}")
plt.close()

# Plotting Density
OUTPUT_FIG_RHO = os.path.join(FIG_DIR, 'density_profile.png')
rho_sim = final_data[:, 4]
plt.figure(figsize=(8, 6))
plt.plot(rho_sim, y_sim, 'ms-', markerfacecolor='none', markersize=6, label=f'LBM Simulation (t={int(last_timestep)})')
plt.xlabel('Density (Lattice Units)')
plt.ylabel('$y$ (Lattice Units)')
plt.title('Steady-State Density Profile')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
# Adjust xlim if needed based on expected density variation
plt.xlim(rho0 - 0.01, rho0 + 0.01) # Example zoom around rho0

plt.savefig(OUTPUT_FIG_RHO)
print(f"Density profile plot saved to {OUTPUT_FIG_RHO}")
plt.close() 