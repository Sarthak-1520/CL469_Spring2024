import numpy as np
import matplotlib.pyplot as plt
import os

# --- Simulation Parameters (should match parameters.dat used for the run) ---
# Need H, g, rho0, Nx to calculate non-dimensional force scale
# TODO: Read these from parameters.dat or results_summary.txt for robustness
H = 11      # Channel height in lattice units (Ny)
Nx = H      # Channel length in lattice units (Assuming Lx=H)
g = 9.203e-6 # Gravity (lattice units)
rho0 = 1.0  # Reference density
# --------------------------------------------------------------------------

DATA_DIR = 'data/' # Use path relative to assignment_two root
FIG_DIR = 'figures/'
FORCE_FILE = os.path.join(DATA_DIR, 'forces.dat')
OUTPUT_FIG = os.path.join(FIG_DIR, 'forces_vs_time.png')
OUTPUT_FIG_NONDIM = os.path.join(FIG_DIR, 'forces_nondim_vs_time.png')
OUTPUT_FIG_YFORCE = os.path.join(FIG_DIR, 'forces_y_vs_time.png')
OUTPUT_FIG_XFORCE_COMP = os.path.join(FIG_DIR, 'force_methods_comparison.png')

# Ensure figure directory exists
os.makedirs(FIG_DIR, exist_ok=True)

# Analytical non-dimensional force (Eq. 8)
F_analytical_nondim = 0.5

# Calculate the characteristic force scale (Weight of fluid in channel section)
# W_lattice = rho0 * g * Volume = rho0 * g * (Lx * dx) * (H_width * dy) * (1 * dz)
# Assuming dx=dy=dz=1, Lx=Nx, H_width = Ny-1 = H-1
H_width = H - 1
W_lattice = rho0 * g * Nx * H_width
print(f"Characteristic force scale (Weight in lattice units): {W_lattice:.4e}")

if W_lattice < 1e-15:
    print("Warning: Characteristic force scale is near zero. Non-dimensionalization might be unstable.")
    # Avoid division by zero later
    W_lattice = 1.0

# Load simulation data
try:
    # Load header to map column names to indices
    with open(FORCE_FILE, 'r') as f:
        header_line = f.readline().strip()
        if header_line.startswith('#'):
            header = header_line[1:].strip().split()
            col_map = {name: idx for idx, name in enumerate(header)}
        else:
            # Assume default order if no header
            print("Warning: No header found in forces.dat. Assuming default column order.")
            col_map = {'Timestep': 0, 'F_bot_x_ME': 1, 'F_bot_y_ME': 2, 'F_top_x_ME': 3, 'F_top_y_ME': 4,
                       'F_bot_x_SI': 5, 'F_bot_y_SI': 6, 'F_top_x_SI': 7, 'F_top_y_SI': 8,
                       'F_bot_x_FD': 9, 'F_top_x_FD': 10, 'Dissipation': 11}

    data = np.loadtxt(FORCE_FILE)
except OSError:
    print(f"Error: Could not find or read {FORCE_FILE}")
    print("Please run the simulation first.")
    exit()
except ValueError:
    print(f"Error: Could not parse data in {FORCE_FILE}. Is it complete?")
    exit()
except Exception as e:
     print(f"An unexpected error occurred reading {FORCE_FILE}: {e}")
     exit()

if data.size == 0:
    print(f"Error: {FORCE_FILE} is empty.")
    exit()

# Extract columns based on header mapping
timesteps = data[:, col_map['Timestep']]
F_bot_x_me = data[:, col_map['F_bot_x_ME']]
F_top_x_me = data[:, col_map['F_top_x_ME']]
F_bot_x_si = data[:, col_map['F_bot_x_SI']]
F_top_x_si = data[:, col_map['F_top_x_SI']]
F_bot_x_fd = data[:, col_map['F_bot_x_FD']]
F_top_x_fd = data[:, col_map['F_top_x_FD']]

# Extract Y forces if available
F_bot_y_me = data[:, col_map.get('F_bot_y_ME', -1)] if 'F_bot_y_ME' in col_map else np.zeros_like(timesteps)
F_top_y_me = data[:, col_map.get('F_top_y_ME', -1)] if 'F_top_y_ME' in col_map else np.zeros_like(timesteps)
F_bot_y_si = data[:, col_map.get('F_bot_y_SI', -1)] if 'F_bot_y_SI' in col_map else np.zeros_like(timesteps)
F_top_y_si = data[:, col_map.get('F_top_y_SI', -1)] if 'F_top_y_SI' in col_map else np.zeros_like(timesteps)

# --- Plot Raw Forces (Momentum Exchange) vs Time ---
plt.figure(figsize=(10, 6))
plt.plot(timesteps, F_bot_x_me, 'b-', label='Bottom Wall Force (ME, x)')
plt.plot(timesteps, F_top_x_me, 'r--', label='Top Wall Force (ME, x)')

plt.xlabel('Timestep')
plt.ylabel('Force (Lattice Units)')
plt.title('Wall Forces (Momentum Exchange) vs. Time')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

plt.savefig(OUTPUT_FIG)
print(f"Raw force plot saved to {OUTPUT_FIG}")
plt.close()

# --- Plot Non-Dimensional Forces vs Time ---
F_bot_x_me_nd = F_bot_x_me / W_lattice
F_top_x_me_nd = F_top_x_me / W_lattice
F_bot_x_si_nd = F_bot_x_si / W_lattice
F_top_x_si_nd = F_top_x_si / W_lattice
F_bot_x_fd_nd = F_bot_x_fd / W_lattice
F_top_x_fd_nd = F_top_x_fd / W_lattice

plt.figure(figsize=(10, 6))
plt.plot(timesteps, F_bot_x_me_nd, label='Bottom Wall (ME)', color='blue', linestyle='-')
plt.plot(timesteps, F_top_x_me_nd, label='Top Wall (ME)', color='red', linestyle='-')

plt.plot(timesteps, F_bot_x_si_nd, label='Bottom Wall (SI)', color='cyan', linestyle='--')
plt.plot(timesteps, F_top_x_si_nd, label='Top Wall (SI)', color='magenta', linestyle='--')

plt.plot(timesteps, F_bot_x_fd_nd, label='Bottom Wall (FD)', color='green', linestyle=':')
plt.plot(timesteps, F_top_x_fd_nd, label='Top Wall (FD)', color='orange', linestyle=':')

plt.axhline(F_analytical_nondim, color='black', linestyle='-',
            linewidth=2, label=f'Analytical Value ({F_analytical_nondim})')

plt.xlabel('Timestep')
plt.ylabel('Non-dimensional Force ($F_x / W_{lattice}$)')
plt.title('Non-dimensional Wall Forces vs. Time')
plt.legend(ncol=2)
plt.grid(True, linestyle='--', alpha=0.6)
plt.ylim(bottom=0) # Force should be positive

plt.savefig(OUTPUT_FIG_NONDIM)
print(f"Non-dimensional force plot saved to {OUTPUT_FIG_NONDIM}")
plt.close()

# --- Plot Y-Forces vs Time ---
plt.figure(figsize=(10, 6))
plt.plot(timesteps, F_bot_y_me, label='Bottom Wall (ME, y)', color='blue', linestyle='-')
plt.plot(timesteps, F_top_y_me, label='Top Wall (ME, y)', color='red', linestyle='-')
plt.plot(timesteps, F_bot_y_si, label='Bottom Wall (SI, y)', color='cyan', linestyle='--')
plt.plot(timesteps, F_top_y_si, label='Top Wall (SI, y)', color='magenta', linestyle='--')
plt.xlabel('Timestep')
plt.ylabel('Force Fy (Lattice Units)')
plt.title('Wall Y-Forces vs. Time')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(OUTPUT_FIG_YFORCE)
print(f"Y-force plot saved to {OUTPUT_FIG_YFORCE}")
plt.close()

# --- Plot X-Force Method Comparison vs Time ---
plt.figure(figsize=(10, 6))
# Plot only one wall for clarity, e.g., bottom wall
plt.plot(timesteps, F_bot_x_me, label='Bottom Wall (ME)', color='blue', linestyle='-')
plt.plot(timesteps, F_bot_x_si, label='Bottom Wall (SI)', color='cyan', linestyle='--')
plt.plot(timesteps, F_bot_x_fd, label='Bottom Wall (FD)', color='green', linestyle=':')
plt.xlabel('Timestep')
plt.ylabel('Force Fx (Lattice Units)')
plt.title('Comparison of X-Force Calculation Methods (Bottom Wall)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(OUTPUT_FIG_XFORCE_COMP)
print(f"X-force method comparison plot saved to {OUTPUT_FIG_XFORCE_COMP}")
plt.close()

# Print final non-dimensional forces
print(f"\nFinal Non-dimensional Forces (at t={int(timesteps[-1])}):")
print(f"  Analytical Target: {F_analytical_nondim:.4f}")
print(f"  Bottom ME: {F_bot_x_me_nd[-1]:.4f}")
print(f"  Top ME:    {F_top_x_me_nd[-1]:.4f}")
print(f"  Bottom SI: {F_bot_x_si_nd[-1]:.4f}")
print(f"  Top SI:    {F_top_x_si_nd[-1]:.4f}")
print(f"  Bottom FD: {F_bot_x_fd_nd[-1]:.4f}")
print(f"  Top FD:    {F_top_x_fd_nd[-1]:.4f}")

# Print final raw forces for comparison
F_analytical_lattice = F_analytical_nondim * W_lattice
print(f"\nFinal Raw Forces (Lattice Units, at t={int(timesteps[-1])}):")
print(f"  Analytical Target: {F_analytical_lattice:.4e}")
print(f"  Bottom ME: {F_bot_x_me[-1]:.4e}")
print(f"  Top ME:    {F_top_x_me[-1]:.4e}")
print(f"  Bottom SI: {F_bot_x_si[-1]:.4e}")
print(f"  Top SI:    {F_top_x_si[-1]:.4e}")
print(f"  Bottom FD: {F_bot_x_fd[-1]:.4e}")
print(f"  Top FD:    {F_top_x_fd[-1]:.4e}") 