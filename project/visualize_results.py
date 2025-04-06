#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import glob
import os
import re
from pathlib import Path

# --- Configuration ---
RESULTS_DIR = Path("./results")
CONVERGENCE_FILE = RESULTS_DIR / "convergence.dat"
CONVERGENCE_PLOT_FILE = RESULTS_DIR / "convergence_py.png"
DAT_PATTERN = str(RESULTS_DIR / "lbm_output_*.dat")
ANIMATION_FILE = RESULTS_DIR / "simulation_animation.gif" # Or .mp4 if ffmpeg is installed
SINGLE_PLOT_FILE_TPL = RESULTS_DIR / "field_plot_{timestep}_py.png"

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


# --- Main Execution ---
if __name__ == "__main__":
    print("--- Starting Python Visualization ---")

    # 1. Plot Convergence
    plot_convergence(CONVERGENCE_FILE, CONVERGENCE_PLOT_FILE)

    # 2. Find DAT files
    dat_files = sorted(glob.glob(DAT_PATTERN), key=lambda f: int(re.search(r'_(\d+)\.dat$', f).group(1))) # Sort numerically

    if not dat_files:
        print("No DAT files found to visualize.")
    else:
        # 3. Plot the latest field snapshot
        latest_dat_file = dat_files[-1]
        match = re.search(r'_(\d+)\.dat$', latest_dat_file)
        timestep = match.group(1) if match else "latest"
        field_plot_filename = SINGLE_PLOT_FILE_TPL.with_name(f"field_plot_{timestep}_py.png")
        plot_fields_mpl(latest_dat_file, field_plot_filename, timestep)

        # 4. Create and save animation
        create_animation(dat_files, ANIMATION_FILE)


    print("--- Python Visualization Finished ---")