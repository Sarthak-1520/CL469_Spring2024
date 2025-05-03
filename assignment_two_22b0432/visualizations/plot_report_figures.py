import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import pandas as pd # Using pandas for easier data handling with headers

# --- Constants ---
DATA_DIR_BASE = 'data' # Base directory relative to assignment_two root
FIG_DIR = 'figures/'
ANALYTICAL_NONDIM_FORCE = 0.5

# Ensure figure directory exists
os.makedirs(FIG_DIR, exist_ok=True)

def read_params_from_summary(summary_file):
    """Reads key parameters from the results_summary.txt file."""
    params = {}
    try:
        with open(summary_file, 'r') as f:
            for line in f:
                if line.startswith('#') or ':' not in line:
                    continue
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()
                try:
                    if key == 'H' or key == 'Nx' or key == 'Ny':
                        params[key] = int(value)
                    elif key == 'Collision Operator':
                         params[key] = value
                    else:
                         params[key] = float(value)
                except ValueError:
                    params[key] = value # Store as string if not convertible
    except FileNotFoundError:
        print(f"Warning: Could not find summary file {summary_file}")
    return params

def calculate_nondim_factor(params):
    """ Calculates the factor W_lattice to non-dimensionalize forces."""
    H = params.get('H')
    Nx = params.get('Nx')
    g = params.get('g')
    rho0 = params.get('rho0', 1.0) # Assume rho0=1 if not in summary

    if H is None or Nx is None or g is None:
        print("Warning: Missing parameters (H, Nx, g) in summary file for non-dimensionalization.")
        return None

    H_width = H - 1
    W_lattice = rho0 * g * Nx * H_width
    if W_lattice < 1e-15:
        print("Warning: Characteristic force scale W_lattice is near zero.")
        return None
    return W_lattice

def calculate_force_error(sim_force_nd):
    """Calculates relative percentage error from analytical value 0.5."""
    return (sim_force_nd - ANALYTICAL_NONDIM_FORCE) / ANALYTICAL_NONDIM_FORCE * 100.0

def plot_force_error_vs_resolution(H_values, data_dirs, fig_path):
    """Generates plots for report section 5."""
    results = {
        'H': [],
        'ME_bot_err': [], 'ME_top_err': [],
        'SI_bot_err': [], 'SI_top_err': [],
        'FD_bot_err': [], 'FD_top_err': []
    }

    print("\n--- Generating Force Error vs Resolution Plot ---")
    for H, data_dir in zip(H_values, data_dirs):
        print(f"Processing H={H} from directory {data_dir}...")
        summary_file = os.path.join(data_dir, 'results_summary.txt')
        force_file = os.path.join(data_dir, 'forces.dat')

        params = read_params_from_summary(summary_file)
        if not params or params.get('H') != H:
            print(f"  Skipping H={H}: Summary file missing or H mismatch.")
            continue

        W_lattice = calculate_nondim_factor(params)
        if W_lattice is None:
            print(f"  Skipping H={H}: Could not calculate W_lattice.")
            continue

        try:
            # df = pd.read_csv(force_file, delim_whitespace=True, comment='#')
            # Read header separately and then data
            with open(force_file, 'r') as f:
                header_line = f.readline().strip()
            if header_line.startswith('#'):
                header_cols = header_line[1:].strip().split()
            else:
                # Fallback if header is missing (should not happen based on C++ code)
                header_cols = None
            
            df = pd.read_csv(force_file, sep='\s+', comment='#', names=header_cols, header=0)

            if df.empty:
                 print(f"  Skipping H={H}: Force file is empty.")
                 continue
            final_forces = df.iloc[-1]
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f"  Skipping H={H}: Force file missing or empty.")
            continue
        except Exception as e:
            print(f"  Skipping H={H}: Error reading force file: {e}")
            continue

        results['H'].append(H)
        results['ME_bot_err'].append(calculate_force_error(final_forces['F_bot_x_ME'] / W_lattice))
        results['ME_top_err'].append(calculate_force_error(final_forces['F_top_x_ME'] / W_lattice))
        results['SI_bot_err'].append(calculate_force_error(final_forces['F_bot_x_SI'] / W_lattice))
        results['SI_top_err'].append(calculate_force_error(final_forces['F_top_x_SI'] / W_lattice))
        results['FD_bot_err'].append(calculate_force_error(final_forces['F_bot_x_FD'] / W_lattice))
        results['FD_top_err'].append(calculate_force_error(final_forces['F_top_x_FD'] / W_lattice))

    if not results['H']:
        print("Error: No valid data found for any resolution to plot force errors.")
        return

    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    H_plot = results['H']

    # Bottom Wall Errors
    axes[0].plot(H_plot, results['ME_bot_err'], 'bo-', label='Momentum Exch.')
    axes[0].plot(H_plot, results['SI_bot_err'], 'gs--', label='Stress Integr.')
    axes[0].plot(H_plot, results['FD_bot_err'], 'r^:', label='Finite Diff.')
    axes[0].set_xlabel('Channel Width H (lattice units)')
    axes[0].set_ylabel('Relative Error in $F_{x, nd}$ (%)')
    axes[0].set_title('Bottom Wall Force Error')
    axes[0].grid(True, linestyle='--', alpha=0.6)
    axes[0].legend()

    # Top Wall Errors
    axes[1].plot(H_plot, results['ME_top_err'], 'bo-', label='Momentum Exch.')
    axes[1].plot(H_plot, results['SI_top_err'], 'gs--', label='Stress Integr.')
    axes[1].plot(H_plot, results['FD_top_err'], 'r^:', label='Finite Diff.')
    axes[1].set_xlabel('Channel Width H (lattice units)')
    # axes[1].set_ylabel('Relative Error (%)') # Shared Y
    axes[1].set_title('Top Wall Force Error')
    axes[1].grid(True, linestyle='--', alpha=0.6)
    axes[1].legend()

    plt.suptitle('Force Calculation Method Accuracy vs. Resolution')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout for suptitle
    plt.savefig(fig_path)
    print(f"Force error plot saved to {fig_path}")
    plt.close()

def plot_slip_velocity(tau_values, data_dirs_bgk, data_dirs_trt, fig_path):
    """Generates plots for report section 6 (slip velocity)."""
    print("\n--- Generating Slip Velocity Plot ---")
    # Create figure with two subplots (bottom and top walls)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    colors = plt.cm.viridis(np.linspace(0, 0.8, len(tau_values)))

    # --- Process BGK --- 
    for i, tau in enumerate(tau_values):
        data_dir = data_dirs_bgk.get(tau)
        if not data_dir:
            print(f"  Skipping BGK tau={tau}: No data directory specified.")
            continue
        print(f"Processing BGK tau={tau} from directory {data_dir}...")
        vel_file = os.path.join(data_dir, 'velocity_profile.dat')
        try:
            # df = pd.read_csv(vel_file, delim_whitespace=True, comment='#')
            # Explicitly handle header for velocity file too
            with open(vel_file, 'r') as f:
                header_line = f.readline().strip()
            if header_line.startswith('#'):
                header_cols = header_line[1:].strip().split()
            else:
                header_cols = ['Timestep', 'Y', 'UX', 'UY', 'RHO'] # Default assumed order
            
            df = pd.read_csv(vel_file, sep='\s+', comment='#', names=header_cols, header=0)
            
            if df.empty:
                print(f"  Skipping BGK tau={tau}: Velocity file empty.")
                continue
            final_vel = df[df['Timestep'] == df['Timestep'].iloc[-1]].copy()
            H_sim = final_vel['Y'].max() # Get actual H from data
            final_vel['Y_norm'] = final_vel['Y'] / H_sim # Normalize Y

            # Bottom Wall (y=0)
            axes[0].plot(final_vel['UX'].iloc[0:5], final_vel['Y'].iloc[0:5],
                     marker='o', linestyle='-', color=colors[i],
                     label=f'BGK $\\tau={tau}$')
            # Top Wall (y=H)
            axes[1].plot(final_vel['UX'].iloc[-5:], final_vel['Y'].iloc[-5:],
                     marker='o', linestyle='-', color=colors[i],
                     label=f'BGK $\\tau={tau}$')
        except Exception as e:
            print(f"  Skipping BGK tau={tau}: Error reading/plotting velocity file: {e}")

    # --- Process TRT --- 
    for i, tau_plus in enumerate(tau_values): # Using tau_values as tau_plus for simplicity
        data_dir = data_dirs_trt.get(tau_plus)
        if not data_dir:
             print(f"  Skipping TRT tau+={tau_plus}: No data directory specified.")
             continue
        print(f"Processing TRT tau+={tau_plus} from directory {data_dir}...")
        vel_file = os.path.join(data_dir, 'velocity_profile.dat')
        try:
            df = pd.read_csv(vel_file, delim_whitespace=True, comment='#')
            if df.empty:
                 print(f"  Skipping TRT tau+={tau_plus}: Velocity file empty.")
                 continue
            final_vel = df[df['Timestep'] == df['Timestep'].iloc[-1]].copy()
            H_sim = final_vel['Y'].max() # Get actual H from data
            final_vel['Y_norm'] = final_vel['Y'] / H_sim # Normalize Y
            
            # Bottom Wall (y=0)
            axes[0].plot(final_vel['UX'].iloc[0:5], final_vel['Y'].iloc[0:5],
                     marker='s', linestyle='--', color=colors[i],
                     label=f'TRT $\\tau^+={tau_plus}$') # Assumes tau- is set correctly in run
            # Top Wall (y=H)
            axes[1].plot(final_vel['UX'].iloc[-5:], final_vel['Y'].iloc[-5:],
                     marker='s', linestyle='--', color=colors[i],
                     label=f'TRT $\\tau^+={tau_plus}$')
        except Exception as e:
            print(f"  Skipping TRT tau+={tau_plus}: Error reading/plotting velocity file: {e}")

    # Formatting for Bottom Wall subplot
    axes[0].axvline(0, color='k', linestyle=':', label='No Slip (u=0)')
    axes[0].set_xlabel('Velocity $u_x$ (lattice units)')
    axes[0].set_ylabel('Y coordinate (lattice units)')
    axes[0].set_title('Near Bottom Wall (y=0)')
    axes[0].legend()
    axes[0].grid(True, linestyle='--', alpha=0.6)
    axes[0].set_ylim(-0.5, 4.5) # Zoom near bottom wall

    # Formatting for Top Wall subplot
    axes[1].axvline(0, color='k', linestyle=':', label='No Slip (u=0)')
    axes[1].set_xlabel('Velocity $u_x$ (lattice units)')
    # axes[1].set_ylabel('Y coordinate (lattice units)') # Shared Y axis
    axes[1].set_title(f'Near Top Wall (y={H_sim:.0f})') # Use H from data
    axes[1].legend()
    axes[1].grid(True, linestyle='--', alpha=0.6)
    axes[1].set_ylim(H_sim - 4.5, H_sim + 0.5) # Zoom near top wall

    plt.suptitle('Near-Wall Velocity Profile Comparison (BGK vs TRT)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(fig_path)
    print(f"Slip velocity plot saved to {fig_path}")
    plt.close()

def plot_dissipation_vs_resolution(H_values, data_dirs, fig_path):
    """Generates plots for report section 7 (viscous dissipation)."""
    # Analytical dissipation Wd,p = 1/12 * rho_p^2 * g_p^2 * Ls * H_p^3 * dzp
    # Need to find analytical dissipation in lattice units and non-dimensionalize
    # Wd_lattice = Sum_V (tau - 0.5)/(rho*cs2*tau^2) * Pi_xy^2
    # Analytical Pi_xy = -(tau - 0.5) * rho * cs2 * (dux/dy)
    # dux/dy = g/(2*nu) * (H_width - 2*y)
    # Pi_xy = -(tau-0.5)*rho*cs2 * g/(2*nu) * (H_width - 2*y)
    # Substitute into Wd formula... complicated.
    # Alternative: Use analytical Wd,p and convert to lattice units, or non-dim both.
    # Let's calculate non-dimensional dissipation Wd* = Wd_sim / (rho0 * um^2 * nu)? No.
    # Compare Wd_sim (from Eq 37) with analytical Wd_lattice (from Eq 39 conversion)
    # Wd,p = (1/12) * rho_p^2 * g_p^2 * Ls * H_p^3 * dzp / mu_p  (Eq 39 seems off, mu missing?)
    # Let's use Wd = Integral( mu * (dux/dy)^2 dV )
    # dux/dy = g/(2*nu) * (H_width - 2y)
    # (dux/dy)^2 = (g/(2*nu))^2 * (H_width - 2y)^2
    # Integral_y=0^Hw mu * (g/(2*nu))^2 * (H_width - 2y)^2 dy * (Lx*dz)
    # = mu * (g/(2*nu))^2 * Lx * dz * [ Hw^2*y - 2*Hw*y^2 + 4/3*y^3 ]_0^Hw
    # = mu * (g/(2*nu))^2 * Lx * dz * (Hw^3 - 2*Hw^3 + 4/3*Hw^3)
    # = mu * (g/(2*nu))^2 * Lx * dz * (1/3 * Hw^3)
    # = (rho*nu) * g^2 / (4*nu^2) * Lx * dz * (1/3 * Hw^3)
    # = (rho * g^2 * Lx * dz * Hw^3) / (12 * nu)
    results = {
        'H': [],
        'Wd_sim': [],
        'Wd_analytical': []
    }
    print("\n--- Generating Viscous Dissipation vs Resolution Plot ---")
    for H, data_dir in zip(H_values, data_dirs):
        print(f"Processing H={H} from directory {data_dir}...")
        summary_file = os.path.join(data_dir, 'results_summary.txt')
        force_file = os.path.join(data_dir, 'forces.dat')

        params = read_params_from_summary(summary_file)
        if not params or params.get('H') != H:
            print(f"  Skipping H={H}: Summary file missing or H mismatch.")
            continue

        try:
            # df = pd.read_csv(force_file, delim_whitespace=True, comment='#')
            # Explicitly handle header for dissipation file (which is forces.dat)
            with open(force_file, 'r') as f:
                header_line = f.readline().strip()
            if header_line.startswith('#'):
                header_cols = header_line[1:].strip().split()
            else:
                header_cols = None # Fallback, though header should exist
            
            df = pd.read_csv(force_file, sep='\s+', comment='#', names=header_cols, header=0)

            if df.empty:
                 print(f"  Skipping H={H}: Dissipation file empty.")
                 continue
            final_dissipation = df['Dissipation'].iloc[-1]
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f"  Skipping H={H}: Dissipation file missing or empty.")
            continue
        except Exception as e:
            print(f"  Skipping H={H}: Error reading dissipation file: {e}")
            continue

        # Calculate analytical dissipation in lattice units
        rho0 = params.get('rho0', 1.0)
        g_lattice = params.get('g')
        nu_lattice = params.get('nu')
        Nx_lattice = params.get('Nx')
        Hw_lattice = H - 1

        if g_lattice is None or nu_lattice is None or Nx_lattice is None or nu_lattice < 1e-12:
            print(f"  Skipping H={H}: Missing parameters for analytical dissipation.")
            continue

        Wd_analytical_lattice = (rho0 * g_lattice**2 * Nx_lattice * Hw_lattice**3) / (12.0 * nu_lattice)

        results['H'].append(H)
        results['Wd_sim'].append(final_dissipation)
        results['Wd_analytical'].append(Wd_analytical_lattice)

    if not results['H']:
        print("Error: No valid data found for any resolution to plot dissipation.")
        return

    # Plotting Absolute Values
    plt.figure(figsize=(8, 6))
    H_plot = results['H']
    Wd_sim_plot = results['Wd_sim']
    Wd_analytical_plot = results['Wd_analytical']

    plt.plot(H_plot, Wd_sim_plot, 'bo-', label='Simulation (Eq. 37)')
    plt.plot(H_plot, Wd_analytical_plot, 'r--', label='Analytical (Integral)')

    plt.xlabel('Channel Width H (lattice units)')
    plt.ylabel('Total Viscous Dissipation $W_d$ (lattice units)')
    plt.title('Viscous Dissipation vs. Resolution')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(fig_path)
    print(f"Dissipation plot saved to {fig_path}")
    plt.close()

    # Plotting Relative Error
    OUTPUT_FIG_DISS_ERR = os.path.join(FIG_DIR, 'dissipation_error_vs_resolution.png')
    # Avoid division by zero if analytical is zero (shouldn't happen here)
    Wd_analytical_plot_np = np.array(Wd_analytical_plot)
    valid_idx = Wd_analytical_plot_np != 0
    if np.any(valid_idx):
      error_percent = np.full_like(Wd_analytical_plot_np, np.nan)
      error_percent[valid_idx] = (np.array(Wd_sim_plot)[valid_idx] - Wd_analytical_plot_np[valid_idx]) / Wd_analytical_plot_np[valid_idx] * 100.0
      
      plt.figure(figsize=(8, 6))
      plt.plot(np.array(H_plot)[valid_idx], error_percent[valid_idx], 'gd:', label='Relative Error')
      plt.xlabel('Channel Width H (lattice units)')
      plt.ylabel('Relative Error in $W_d$ (%)')
      plt.title('Viscous Dissipation Error vs. Resolution')
      plt.axhline(0, color='k', linestyle='--', alpha=0.5)
      plt.legend()
      plt.grid(True, linestyle='--', alpha=0.6)
      plt.savefig(OUTPUT_FIG_DISS_ERR)
      print(f"Dissipation error plot saved to {OUTPUT_FIG_DISS_ERR}")
      plt.close()
    else:
      print("Skipping dissipation error plot: Analytical dissipation is zero or no valid data.")

# --- Main Execution --- (Example Usage)
if __name__ == "__main__":
    # These lists define which simulation results to use for the comparison plots.
    # You need to run the simulation multiple times with different parameters
    # and save the results in appropriately named directories.

    # --- Example for Force Error vs Resolution (Section 5) ---
    H_resolution_study = [5, 10, 20, 40] # Example H values
    # Assume data is saved in folders like data_H5/, data_H10/, etc.
    data_dirs_resolution = [os.path.join(DATA_DIR_BASE, f"data_H{h}") for h in H_resolution_study]
    plot_force_error_vs_resolution(H_resolution_study, data_dirs_resolution,
                                   os.path.join(FIG_DIR, 'force_error_vs_resolution.png'))

    # --- Example for Slip Velocity (Section 6) ---
    tau_slip_study = [0.8, 2.0, 5.0] # Example tau values (tau = tau+ for TRT here)
    # Assume data saved in data_BGK_tau0.8/, data_TRT_tau2.0/, etc.
    data_dirs_bgk_slip = {tau: os.path.join(DATA_DIR_BASE, f"data_BGK_tau{tau:.1f}") for tau in tau_slip_study}
    data_dirs_trt_slip = {tau: os.path.join(DATA_DIR_BASE, f"data_TRT_tau{tau:.1f}") for tau in tau_slip_study}
    plot_slip_velocity(tau_slip_study, data_dirs_bgk_slip, data_dirs_trt_slip,
                       os.path.join(FIG_DIR, 'slip_velocity_comparison.png'))

    # --- Example for Dissipation vs Resolution (Section 7) ---
    # Reuse resolution study directories
    plot_dissipation_vs_resolution(H_resolution_study, data_dirs_resolution,
                                   os.path.join(FIG_DIR, 'dissipation_vs_resolution.png'))

    print("\nReport figure generation attempted.")
    print("Please ensure simulation data exists in the expected directories:")
    print(f"  Resolution study dirs: {data_dirs_resolution}")
    print(f"  BGK slip study dirs: {list(data_dirs_bgk_slip.values())}")
    print(f"  TRT slip study dirs: {list(data_dirs_trt_slip.values())}") 