import numpy as np
import matplotlib.pyplot as plt
import os

DATA_DIR = 'data/'
FIG_DIR = 'figures/'
CONVERGENCE_FILE = os.path.join(DATA_DIR, 'convergence.dat')
OUTPUT_FIG = os.path.join(FIG_DIR, 'convergence.png')

# Ensure figure directory exists
os.makedirs(FIG_DIR, exist_ok=True)

# Load simulation data
try:
    data = np.loadtxt(CONVERGENCE_FILE)
except OSError:
    print(f"Error: Could not find or read {CONVERGENCE_FILE}")
    print("Please run the simulation first.")
    exit()
except ValueError:
    print(f"Error: Could not parse data in {CONVERGENCE_FILE}. Is it complete?")
    exit()

if data.size == 0:
    print(f"Error: {CONVERGENCE_FILE} is empty.")
    exit()

# Extract columns
timesteps = data[:, 0]
rel_vel_change = data[:, 1]

# Plotting
plt.figure(figsize=(10, 6))
plt.semilogy(timesteps, rel_vel_change, 'b-') # Use log scale for y-axis

plt.xlabel('Timestep')
plt.ylabel('Relative Velocity Change (L2 Norm)')
plt.title('Simulation Convergence')
plt.grid(True, which='both', linestyle='--', alpha=0.6) # Grid for log scale

plt.savefig(OUTPUT_FIG)
print(f"Convergence plot saved to {OUTPUT_FIG}")
plt.close() 