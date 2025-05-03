#!/usr/bin/env python3
"""
Simple visualization script for Entropic Lattice Boltzmann Method simulation results.
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re
import argparse
from pathlib import Path

def get_paths(test_case, collision_model):
    """Get all the paths based on test case and collision model."""
    results_dir = Path("./results")
    figures_dir = Path("./figures")
    
    case_dir = results_dir / f"{test_case}_{collision_model}"
    convergence_file = case_dir / "convergence.dat"
    dat_pattern = str(case_dir / f"{test_case}_{collision_model}_*.dat")
    
    return {
        'results_dir': results_dir,
        'figures_dir': figures_dir,
        'case_dir': case_dir,
        'convergence_file': convergence_file,
        'dat_pattern': dat_pattern
    }

def read_dat_file(filename):
    """Reads the custom .dat file format."""
    try:
        with open(filename, 'r') as f:
            # Read header line 1 (# NX NY)
            f.readline()
            # Read NX NY values
            nx, ny = map(int, f.readline().split())
            # Read header line 3 (# i j rho u.x u.y is_fluid)
            f.readline()
            # Read data using numpy, skipping comment lines is implicit
            data = np.loadtxt(f)

        # Reshape data based on nx, ny assuming row-major (y varies faster) output from C++
        # Columns: 0:i, 1:j, 2:rho, 3:ux, 4:uy, 5:is_fluid
        # We need rho[i, j], ux[i, j], uy[i, j]
        # Create grids for coordinates
        x_coords = data[:, 0].reshape(ny, nx)
        y_coords = data[:, 1].reshape(ny, nx)
        # Extract and reshape data fields
        rho = data[:, 2].reshape(ny, nx)
        ux = data[:, 3].reshape(ny, nx)
        uy = data[:, 4].reshape(ny, nx)
        is_fluid = data[:, 5].reshape(ny, nx).astype(bool)

        # Apply mask based on is_fluid (set non-fluid values to NaN for plotting)
        rho[~is_fluid] = np.nan
        ux[~is_fluid] = np.nan
        uy[~is_fluid] = np.nan

        return x_coords, y_coords, rho, ux, uy, nx, ny
    except Exception as e:
        print(f"Error reading DAT file {filename}: {e}")
        return None, None, None, None, None, 0, 0

def calculate_vorticity(ux, uy, dx=1.0, dy=1.0):
    """Calculate vorticity from velocity field using central differences."""
    # Initialize vorticity array
    vorticity = np.zeros_like(ux)

    # Use central differences for interior points
    # vorticity = duy/dx - dux/dy
    vorticity[1:-1, 1:-1] = (uy[1:-1, 2:] - uy[1:-1, :-2]) / (2 * dx) - \
                            (ux[2:, 1:-1] - ux[:-2, 1:-1]) / (2 * dy)

    # Forward/backward differences for boundaries
    # Left/right boundaries
    vorticity[1:-1, 0] = (uy[1:-1, 1] - uy[1:-1, 0]) / dx - \
                         (ux[2:, 0] - ux[:-2, 0]) / (2 * dy)
    vorticity[1:-1, -1] = (uy[1:-1, -1] - uy[1:-1, -2]) / dx - \
                          (ux[2:, -1] - ux[:-2, -1]) / (2 * dy)

    # Top/bottom boundaries
    vorticity[0, 1:-1] = (uy[0, 2:] - uy[0, :-2]) / (2 * dx) - \
                         (ux[1, 1:-1] - ux[0, 1:-1]) / dy
    vorticity[-1, 1:-1] = (uy[-1, 2:] - uy[-1, :-2]) / (2 * dx) - \
                          (ux[-1, 1:-1] - ux[-2, 1:-1]) / dy

    # Corners
    vorticity[0, 0] = (uy[0, 1] - uy[0, 0]) / dx - \
                      (ux[1, 0] - ux[0, 0]) / dy
    vorticity[0, -1] = (uy[0, -1] - uy[0, -2]) / dx - \
                       (ux[1, -1] - ux[0, -1]) / dy
    vorticity[-1, 0] = (uy[-1, 1] - uy[-1, 0]) / dx - \
                       (ux[-1, 0] - ux[-2, 0]) / dy
    vorticity[-1, -1] = (uy[-1, -1] - uy[-1, -2]) / dx - \
                        (ux[-1, -1] - ux[-2, -1]) / dy

    return vorticity

def plot_velocity_field(dat_filename, output_file, timestep):
    """Reads a DAT file and plots velocity field."""
    print(f"Plotting velocity field from {dat_filename}...")

    x, y, rho, ux, uy, nx, ny = read_dat_file(dat_filename)
    if x is None: return # Error reading file

    # Calculate velocity magnitude
    vel_mag = np.sqrt(ux**2 + uy**2)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8 * (ny / nx) if nx > 0 else 8))

    # Plot velocity magnitude contours
    cont = ax.contourf(x, y, vel_mag, cmap='viridis', levels=50)
    plt.colorbar(cont, ax=ax, label='Velocity Magnitude')

    # Add streamlines
    ux_stream = np.nan_to_num(ux)
    uy_stream = np.nan_to_num(uy)
    ax.streamplot(x, y, ux_stream, uy_stream,
                 color='white', linewidth=0.8, density=1.5,
                 arrowstyle='->', arrowsize=1.0)

    ax.set_title(f"Velocity Field (Step {timestep})")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Velocity field plot saved to {output_file}")
    plt.close(fig)

def plot_vorticity_field(dat_filename, output_file, timestep):
    """Reads a DAT file and plots vorticity field."""
    print(f"Plotting vorticity field from {dat_filename}...")

    x, y, rho, ux, uy, nx, ny = read_dat_file(dat_filename)
    if x is None: return # Error reading file

    # Calculate vorticity
    vorticity = calculate_vorticity(ux, uy)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8 * (ny / nx) if nx > 0 else 8))

    # Create a symmetric colormap (blue-white-red)
    cmap = plt.cm.get_cmap('RdBu_r')

    # Find symmetric vorticity limits for better visualization
    vmax = np.nanmax(np.abs(vorticity))
    vmin = -vmax

    # Plot vorticity contours
    cont = ax.contourf(x, y, vorticity, cmap=cmap, levels=50,
                      vmin=vmin, vmax=vmax)
    plt.colorbar(cont, ax=ax, label='Vorticity')

    # Add streamlines
    ux_stream = np.nan_to_num(ux)
    uy_stream = np.nan_to_num(uy)
    ax.streamplot(x, y, ux_stream, uy_stream,
                 color='black', linewidth=0.8, density=1.5,
                 arrowstyle='->', arrowsize=1.0)

    ax.set_title(f"Vorticity Field (Step {timestep})")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Vorticity field plot saved to {output_file}")
    plt.close(fig)

def plot_convergence(data_file, output_file):
    """Reads convergence data and plots it."""
    print(f"Plotting convergence from {data_file}...")
    if not data_file.exists():
        print(f"Warning: Convergence data file not found: {data_file}")
        return

    try:
        # Load data, skipping header row starting with '#'
        data = np.loadtxt(data_file, comments='#')
        if data.ndim == 1: # Handle case with only one data point
             data = data.reshape(1, -1)
        if data.shape[0] == 0:
             print("Warning: Convergence data file is empty.")
             return

        steps = data[:, 0]
        max_delta_u = data[:, 1]

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(steps, max_delta_u, marker='.', linestyle='-', markersize=3)
        ax.set_yscale('log')
        ax.set_xlabel("Time Step")
        ax.set_ylabel("Max |Δu|")
        ax.set_title("Convergence History (Max Velocity Change)")
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        print(f"Convergence plot saved to {output_file}")
        plt.close(fig)

    except Exception as e:
        print(f"Error plotting convergence: {e}")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Visualize LBM simulation results.')
    parser.add_argument('--type', choices=['velocity', 'vorticity', 'convergence', 'all'],
                        default='all', help='Type of plot to generate')
    parser.add_argument('--step', type=int, default=1000, help='Specific time step to visualize')
    parser.add_argument('--test', type=str, default='ellipse_flow', help='Test case name')
    parser.add_argument('--model', type=str, default='entropic', help='Collision model (entropic or trt)')
    args = parser.parse_args()

    # Get paths based on test case and model
    paths = get_paths(args.test, args.model)
    
    # Create figures directory if it doesn't exist
    os.makedirs(paths['figures_dir'], exist_ok=True)
    
    # Print configuration
    print(f"Processing test case: {args.test} with model: {args.model}")
    print(f"Data directory: {paths['case_dir']}")
    print(f"Convergence file: {paths['convergence_file']}")
    print(f"Data file pattern: {paths['dat_pattern']}")

    # Find the latest timestep file if step is not specified
    if args.step is None:
        dat_files = sorted(glob.glob(paths['dat_pattern']))
        if not dat_files:
            print(f"No data files found matching pattern: {paths['dat_pattern']}")
            return
        latest_file = dat_files[-1]
        match = re.search(r'_(\d+)\.dat$', latest_file)
        if match:
            args.step = int(match.group(1))
            print(f"Using latest timestep: {args.step}")
        else:
            args.step = 1000
            print(f"Could not determine latest timestep, using default: {args.step}")

    # Construct the specific data file path for the requested timestep
    dat_file = paths['case_dir'] / f"{args.test}_{args.model}_{args.step}.dat"
    
    # Check if the file exists
    if not dat_file.exists():
        print(f"Warning: Data file not found: {dat_file}")
        # Try to find the closest timestep
        dat_files = sorted(glob.glob(paths['dat_pattern']))
        if dat_files:
            print(f"Using the latest available timestep instead.")
            dat_file = Path(dat_files[-1])
            match = re.search(r'_(\d+)\.dat$', str(dat_file))
            if match:
                args.step = int(match.group(1))
                print(f"Using timestep: {args.step}")
        else:
            print(f"No data files found matching pattern: {paths['dat_pattern']}")
            return

    # Generate plots based on the type argument
    if args.type in ['velocity', 'all']:
        velocity_file = paths['figures_dir'] / f"velocity_{args.test}_{args.model}_{args.step}.png"
        plot_velocity_field(dat_file, velocity_file, args.step)

    if args.type in ['vorticity', 'all']:
        vorticity_file = paths['figures_dir'] / f"vorticity_{args.test}_{args.model}_{args.step}.png"
        plot_vorticity_field(dat_file, vorticity_file, args.step)

    if args.type in ['convergence', 'all']:
        convergence_file = paths['figures_dir'] / f"convergence_{args.test}_{args.model}.png"
        plot_convergence(paths['convergence_file'], convergence_file)

if __name__ == "__main__":
    main()
