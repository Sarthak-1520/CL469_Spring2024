#!/usr/bin/env python3
"""
Visualization script for Entropic Lattice Boltzmann Method simulation results.

This script reads the output data files from the simulation and creates visualizations
of the flow field, including velocity magnitude, vorticity, and streamlines.
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import numpy as np
import glob
import os
import re
import sys
import argparse
from pathlib import Path

# --- Configuration ---
RESULTS_DIR = Path("./results")
FIGURES_DIR = Path("./figures")

# These will be updated based on command line arguments
TEST_CASE = "ellipse_flow"
COLLISION_MODEL = "entropic"
CASE_DIR = None
CONVERGENCE_FILE = None
DAT_PATTERN = None
ANIMATION_FILE = None
VELOCITY_ANIMATION_FILE = None
VORTICITY_ANIMATION_FILE = None
CONVERGENCE_PLOT_FILE = None
SINGLE_PLOT_FILE_TPL = None

# Initialize paths based on default values
def update_paths():
    global CASE_DIR, CONVERGENCE_FILE, DAT_PATTERN, ANIMATION_FILE
    global VELOCITY_ANIMATION_FILE, VORTICITY_ANIMATION_FILE
    global CONVERGENCE_PLOT_FILE, SINGLE_PLOT_FILE_TPL

    CASE_DIR = RESULTS_DIR / f"{TEST_CASE}_{COLLISION_MODEL}"
    CONVERGENCE_FILE = CASE_DIR / "convergence.dat"
    DAT_PATTERN = str(CASE_DIR / f"{TEST_CASE}_{COLLISION_MODEL}_*.dat")
    ANIMATION_FILE = FIGURES_DIR / f"simulation_{TEST_CASE}_{COLLISION_MODEL}.mp4"
    VELOCITY_ANIMATION_FILE = FIGURES_DIR / f"velocity_{TEST_CASE}_{COLLISION_MODEL}.mp4"
    VORTICITY_ANIMATION_FILE = FIGURES_DIR / f"vorticity_{TEST_CASE}_{COLLISION_MODEL}.mp4"
    CONVERGENCE_PLOT_FILE = FIGURES_DIR / f"convergence_{TEST_CASE}_{COLLISION_MODEL}.png"
    SINGLE_PLOT_FILE_TPL = str(FIGURES_DIR / "{plot_type}_{test_case}_{model}_{timestep}.png")

# Initialize paths with default values
update_paths()

# --- Plot Convergence History ---
def plot_convergence(data_file, output_file):
    """Reads convergence data and plots it using Matplotlib."""
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
        plt.close(fig) # Close figure to free memory

    except Exception as e:
        print(f"Error plotting convergence: {e}")

# --- Read DAT file ---
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

# --- Calculate Vorticity ---
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


# --- Plot LBM Fields (Density and Velocity) using Matplotlib ---
def plot_fields_mpl(dat_filename, output_file, timestep):
    """Reads a DAT file and plots density and velocity using Matplotlib."""
    print(f"Plotting fields from {dat_filename}...")

    x, y, rho, ux, uy, nx, ny = read_dat_file(dat_filename)
    if x is None: return # Error reading file

    fig, ax = plt.subplots(figsize=(8, 8 * (ny / nx) if nx > 0 else 8)) # Adjust aspect ratio

    # Plot density contours
    cont = ax.contourf(x, y, rho, cmap='viridis', levels=50, zorder=1)
    plt.colorbar(cont, ax=ax, label='Density (rho)')

    # Plot velocity vectors (quiver plot) - sample points to avoid clutter
    skip = max(1, nx // 20) # Use a reasonable skip rate for the snapshot
    ax.quiver(x[::skip, ::skip], y[::skip, ::skip],
              ux[::skip, ::skip], uy[::skip, ::skip],
              color='red', # Keep color consistent
              # scale=20,    # Try removing scale here too for consistency
              width=0.003, zorder=10) # Draw quiver on top

    # --- Add Streamlines ---
    # Streamplot uses the *grid* coordinates and the U, V components defined on that grid.
    # Ensure ux and uy don't contain NaNs where streamlines shouldn't go (walls)
    # streamplot handles masked arrays, but filling NaNs with 0 might be safer
    ux_stream = np.nan_to_num(ux)
    uy_stream = np.nan_to_num(uy)
    strm = ax.streamplot(x, y, ux_stream, uy_stream,
                         color='black',      # Streamline color
                         linewidth=0.8,      # Streamline width
                         density=1.5,        # Controls density of streamlines
                         arrowstyle='->',
                         arrowsize=1.0,
                         zorder=5)           # Draw streamlines above contour but below quiver

    # --- End Streamlines ---

    ax.set_title(f"Density, Velocity, and Streamlines (Step {timestep})") # Updated title
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Field plot saved to {output_file}")
    plt.close(fig) # Close figure to free memory

# --- Create Animation ---
def create_animation(dat_files, output_file):
    """Creates an animation from a list of DAT files with a fixed color scale."""
    if not dat_files:
        print("No DAT files found for animation.")
        return

    print(f"Analyzing data range for animation scale...")
    global_rho_min = float('inf')
    global_rho_max = float('-inf')
    valid_files_for_anim = [] # Store files that are successfully read

    for filename in dat_files:
        x, y, rho, ux, uy, nx, ny = read_dat_file(filename)
        if x is not None:
             valid_files_for_anim.append(filename) # Only use valid files
             # Find min/max excluding NaNs (which represent walls)
             current_min = np.nanmin(rho)
             current_max = np.nanmax(rho)
             if not np.isnan(current_min):
                 global_rho_min = min(global_rho_min, current_min)
             if not np.isnan(current_max):
                 global_rho_max = max(global_rho_max, current_max)
        else:
             print(f"Skipping file {filename} due to read error.")

    if not valid_files_for_anim:
         print("No valid DAT files could be read for animation.")
         return
    if global_rho_min == float('inf') or global_rho_max == float('-inf'):
        print("Warning: Could not determine valid global density range. Using default scale.")
        # Fallback or fixed default range if needed
        global_rho_min, global_rho_max = 0.95, 1.05 # Example fallback

    print(f"Global density range for animation: [{global_rho_min:.4f}, {global_rho_max:.4f}]")
    # Define fixed contour levels based on the global range
    contour_levels = np.linspace(global_rho_min, global_rho_max, 51) # 51 levels -> 50 intervals

    print(f"Creating animation from {len(valid_files_for_anim)} frames...")

    # Read the first valid file again to set up the plot
    first_file = valid_files_for_anim[0]
    x, y, rho_init, ux_init, uy_init, nx, ny = read_dat_file(first_file)
    # We already checked read errors above, but double check nx, ny
    if nx == 0 or ny == 0:
        print("Error: Invalid dimensions from first file.")
        return

    fig, ax = plt.subplots(figsize=(8, 8 * (ny / nx))) # Adjust aspect ratio properly
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # Initial plot elements using fixed levels
    cont = ax.contourf(x, y, rho_init, cmap='viridis', levels=contour_levels, vmin=global_rho_min, vmax=global_rho_max)
    cbar = plt.colorbar(cont, ax=ax, label='Density (rho)')
    skip = max(1, nx // 20)
    quiv = ax.quiver(x[::skip, ::skip], y[::skip, ::skip],
                     ux_init[::skip, ::skip], uy_init[::skip, ::skip],
                     color='white', scale=20, width=0.003)
    match = re.search(r'_(\d+)\.dat$', first_file)
    initial_timestep = match.group(1) if match else "0"
    title = ax.set_title(f"Step {initial_timestep}")


    # Update function for animation frames
    def update(frame_num):
        filename = valid_files_for_anim[frame_num] # Use the list of valid files
        match = re.search(r'_(\d+)\.dat$', filename)
        timestep = match.group(1) if match else str(frame_num) # Timestep from filename
        if frame_num % 10 == 0: # Print progress occasionally
            print(f"Processing frame {frame_num+1}/{len(valid_files_for_anim)} (Step {timestep})...")

        _x, _y, rho, ux, uy, _nx, _ny = read_dat_file(filename)
        if _x is None: return [] # Return empty list on error

        # Clear previous plot elements on the axes
        ax.clear()
        # Re-apply static settings
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        current_title = ax.set_title(f"Step {timestep}") # Update title

        # Redraw contour with fixed levels
        current_cont = ax.contourf(_x, _y, rho, cmap='viridis', levels=contour_levels,
                                   vmin=global_rho_min, vmax=global_rho_max, zorder=1)

        # --- Draw Streamlines instead of Quiver ---
        ux_stream = np.nan_to_num(ux) # Ensure no NaNs are passed
        uy_stream = np.nan_to_num(uy)
        # Note: density for streamplot might need tuning for animation speed/clarity
        current_strm = ax.streamplot(_x, _y, ux_stream, uy_stream,
                                     color='black',      # Streamline color
                                     linewidth=0.8,
                                     density=1.5,        # Adjust density if needed
                                     arrowstyle='->',
                                     arrowsize=1.0,
                                     zorder=5)          # Draw streamlines on top of contour
        # --- End Streamlines ---

        # Return artists: contour collections, streamline lines, title
        # streamplot returns a StreamplotSet; we need its lines attribute for blitting (if used)
        # It's generally safer for non-blit animations to just return the core objects
        # Let's return the contour collections, the StreamplotSet's lines, and the title text object
        return list(current_cont.collections) + [current_strm.lines, current_title]


    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=len(valid_files_for_anim),
                                  interval=100, blit=False) # blit=False is safer/easier

    # Save the animation
    try:
        # Try saving as GIF first (requires PillowWriter, often installed with matplotlib)
        writer = animation.PillowWriter(fps=15)
        ani.save(output_file, writer=writer)
        print(f"Animation saved to {output_file}")
    except Exception as e1:
        print(f"Could not save animation as GIF: {e1}")
        print("Try installing Pillow: pip install Pillow")
        # Try saving as MP4 (requires ffmpeg)
        try:
             mp4_output_file = Path(output_file).with_suffix('.mp4')
             writer = animation.FFMpegWriter(fps=15)
             ani.save(str(mp4_output_file), writer=writer)
             print(f"Animation saved to {mp4_output_file}")
        except Exception as e2:
             print(f"Could not save animation as MP4: {e2}")
             print("Ensure ffmpeg is installed and in your PATH.")

    plt.close(fig) # Close figure


# --- Plot Vorticity Field ---
def plot_vorticity_field(dat_filename, output_file, timestep):
    """Reads a DAT file and plots vorticity field using Matplotlib."""
    print(f"Plotting vorticity from {dat_filename}...")

    x, y, rho, ux, uy, nx, ny = read_dat_file(dat_filename)
    if x is None: return # Error reading file

    # Calculate vorticity
    vorticity = calculate_vorticity(ux, uy)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8 * (ny / nx) if nx > 0 else 8)) # Adjust aspect ratio

    # Create a symmetric colormap (blue-white-red)
    cmap = plt.cm.get_cmap('RdBu_r')

    # Find symmetric vorticity limits for better visualization
    vmax = np.nanmax(np.abs(vorticity))
    vmin = -vmax

    # Plot vorticity contours
    cont = ax.contourf(x, y, vorticity, cmap=cmap, levels=50,
                      vmin=vmin, vmax=vmax, zorder=1)
    plt.colorbar(cont, ax=ax, label='Vorticity')

    # Add streamlines
    ux_stream = np.nan_to_num(ux)
    uy_stream = np.nan_to_num(uy)
    ax.streamplot(x, y, ux_stream, uy_stream,
                 color='black',      # Streamline color
                 linewidth=0.8,      # Streamline width
                 density=1.5,        # Controls density of streamlines
                 arrowstyle='->',
                 arrowsize=1.0,
                 zorder=5)           # Draw streamlines above contour

    ax.set_title(f"Vorticity Field (Step {timestep})")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Vorticity plot saved to {output_file}")
    plt.close(fig) # Close figure to free memory

# --- Create Vorticity Animation ---
def create_vorticity_animation(dat_files, output_file):
    """Creates a vorticity animation from a list of DAT files."""
    if not dat_files:
        print("No DAT files found for animation.")
        return

    print(f"Analyzing data range for vorticity animation scale...")
    global_vort_max = 0.0
    valid_files_for_anim = [] # Store files that are successfully read

    for filename in dat_files:
        x, y, rho, ux, uy, nx, ny = read_dat_file(filename)
        if x is not None:
            valid_files_for_anim.append(filename) # Only use valid files
            # Calculate vorticity
            vorticity = calculate_vorticity(ux, uy)
            # Find max vorticity magnitude (for symmetric colormap)
            current_max = np.nanmax(np.abs(vorticity))
            if not np.isnan(current_max):
                global_vort_max = max(global_vort_max, current_max)
        else:
            print(f"Skipping file {filename} due to read error.")

    if not valid_files_for_anim:
        print("No valid DAT files could be read for animation.")
        return
    if global_vort_max == 0.0:
        print("Warning: Could not determine valid vorticity range. Using default scale.")
        global_vort_max = 0.1 # Example fallback

    print(f"Global vorticity range for animation: [-{global_vort_max:.4f}, {global_vort_max:.4f}]")

    # Create figure
    first_file = valid_files_for_anim[0]
    x, y, rho, ux, uy, nx, ny = read_dat_file(first_file)
    if nx == 0 or ny == 0:
        print("Error: Invalid dimensions from first file.")
        return

    fig, ax = plt.subplots(figsize=(8, 8 * (ny / nx))) # Adjust aspect ratio
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # Calculate initial vorticity
    vorticity_init = calculate_vorticity(ux, uy)

    # Create a symmetric colormap (blue-white-red)
    cmap = plt.cm.get_cmap('RdBu_r')

    # Initial plot
    cont = ax.contourf(x, y, vorticity_init, cmap=cmap, levels=50,
                      vmin=-global_vort_max, vmax=global_vort_max)
    cbar = plt.colorbar(cont, ax=ax, label='Vorticity')

    # Add initial streamlines
    ux_stream = np.nan_to_num(ux)
    uy_stream = np.nan_to_num(uy)
    strm = ax.streamplot(x, y, ux_stream, uy_stream,
                        color='black', linewidth=0.8, density=1.5,
                        arrowstyle='->', arrowsize=1.0)

    match = re.search(r'_(\d+)\.dat$', first_file)
    initial_timestep = match.group(1) if match else "0"
    title = ax.set_title(f"Vorticity (Step {initial_timestep})")

    # Update function for animation frames
    def update(frame_num):
        filename = valid_files_for_anim[frame_num]
        match = re.search(r'_(\d+)\.dat$', filename)
        timestep = match.group(1) if match else str(frame_num)
        if frame_num % 10 == 0:
            print(f"Processing vorticity frame {frame_num+1}/{len(valid_files_for_anim)} (Step {timestep})...")

        _x, _y, _rho, ux, uy, _nx, _ny = read_dat_file(filename)
        if _x is None: return []

        # Calculate vorticity
        vorticity = calculate_vorticity(ux, uy)

        # Clear previous plot elements
        ax.clear()
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        current_title = ax.set_title(f"Vorticity (Step {timestep})")

        # Redraw contour with fixed scale
        current_cont = ax.contourf(_x, _y, vorticity, cmap=cmap, levels=50,
                                 vmin=-global_vort_max, vmax=global_vort_max)

        # Add streamlines
        ux_stream = np.nan_to_num(ux)
        uy_stream = np.nan_to_num(uy)
        current_strm = ax.streamplot(_x, _y, ux_stream, uy_stream,
                                   color='black', linewidth=0.8, density=1.5,
                                   arrowstyle='->', arrowsize=1.0)

        return list(current_cont.collections) + [current_strm.lines, current_title]

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(valid_files_for_anim),
                                interval=100, blit=False)

    # Save animation
    try:
        writer = animation.FFMpegWriter(fps=15)
        ani.save(output_file, writer=writer)
        print(f"Vorticity animation saved to {output_file}")
    except Exception as e:
        print(f"Could not save vorticity animation: {e}")
        print("Ensure ffmpeg is installed and in your PATH.")

    plt.close(fig)

# --- Main Execution ---
if __name__ == "__main__":
    print("--- Starting Python Visualization ---")

    # Create figures directory if it doesn't exist
    os.makedirs(FIGURES_DIR, exist_ok=True)

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Visualize LBM simulation results.')
    parser.add_argument('--type', choices=['velocity', 'vorticity', 'streamlines', 'all'],
                        default='all', help='Type of plot to generate')
    parser.add_argument('--animate', action='store_true', help='Create animation')
    parser.add_argument('--convergence', action='store_true', help='Plot convergence history')
    parser.add_argument('--step', type=int, help='Specific time step to visualize')
    parser.add_argument('--test', type=str, default='ellipse_flow', help='Test case name (e.g., ellipse_flow, lid_driven_cavity)')
    parser.add_argument('--model', type=str, default='entropic', help='Collision model (entropic or trt)')
    args = parser.parse_args()

    # Update the global variables with the test case and model
    global TEST_CASE, COLLISION_MODEL

    TEST_CASE = args.test
    COLLISION_MODEL = args.model

    # Update paths based on test case and model
    update_paths()

    # Print configuration
    print(f"Processing test case: {TEST_CASE} with model: {COLLISION_MODEL}")
    print(f"Data directory: {CASE_DIR}")
    print(f"Convergence file: {CONVERGENCE_FILE}")
    print(f"Data file pattern: {DAT_PATTERN}")

    # 1. Plot Convergence
    if args.convergence or args.type == 'all':
        plot_convergence(CONVERGENCE_FILE, CONVERGENCE_PLOT_FILE)

    # 2. Find DAT files
    dat_files = sorted(glob.glob(DAT_PATTERN), key=lambda f: int(re.search(r'_(\d+)\.dat$', f).group(1))) # Sort numerically

    if not dat_files:
        print("No DAT files found to visualize.")
    else:
        # Get the latest time step
        latest_dat_file = dat_files[-1]
        match = re.search(r'_(\d+)\.dat$', latest_dat_file)
        latest_timestep = match.group(1) if match else "latest"

        # Determine which time step to visualize
        timestep = args.step if args.step is not None else latest_timestep
        target_file = f"{RESULTS_DIR}/lbm_output_{timestep}.dat"

        if not os.path.exists(target_file):
            print(f"Warning: File for time step {timestep} not found. Using latest time step {latest_timestep}.")
            target_file = latest_dat_file
            timestep = latest_timestep

        # 3. Plot the field snapshots
        if args.type in ['velocity', 'all']:
            field_plot_filename = SINGLE_PLOT_FILE_TPL.format(
                plot_type="velocity",
                test_case=TEST_CASE,
                model=COLLISION_MODEL,
                timestep=timestep
            )
            plot_fields_mpl(target_file, field_plot_filename, timestep)

        if args.type in ['vorticity', 'all']:
            vorticity_plot_filename = SINGLE_PLOT_FILE_TPL.format(
                plot_type="vorticity",
                test_case=TEST_CASE,
                model=COLLISION_MODEL,
                timestep=timestep
            )
            plot_vorticity_field(target_file, vorticity_plot_filename, timestep)

        # 4. Create and save animations
        if args.animate:
            if args.type in ['velocity', 'all']:
                create_animation(dat_files, ANIMATION_FILE)

            if args.type in ['vorticity', 'all']:
                create_vorticity_animation(dat_files, VORTICITY_ANIMATION_FILE)

    print("--- Python Visualization Finished ---")