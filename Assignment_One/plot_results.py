#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys

def plot_velocity_profile(filename, output_file=None):
    """Plot velocity profile with analytical comparison"""
    print(f"Reading velocity profile from {filename}")
    try:
        data = np.loadtxt(filename)
        print(f"Profile data shape: {data.shape}")
        
        y = data[:,0]
        u_sim = data[:,1]
        u_analytical = data[:,2]
        
        print(f"Y range: {y.min()} to {y.max()}")
        print(f"Simulation velocity range: {u_sim.min()} to {u_sim.max()}")
        print(f"Analytical velocity range: {u_analytical.min()} to {u_analytical.max()}")
        
        plt.figure(figsize=(10, 6))
        plt.plot(u_sim, y, 'o-', label='Simulation')
        plt.plot(u_analytical, y, '--', label='Analytical')
        plt.xlabel('Velocity u$_x$')
        plt.ylabel('Channel height y')
        plt.title('Poiseuille Flow Velocity Profile')
        plt.legend()
        plt.grid(True)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved velocity profile plot to {output_file}")
        else:
            plt.show()
        plt.close()
    except Exception as e:
        print(f"Error in plot_velocity_profile: {e}")
        import traceback
        traceback.print_exc()

def plot_velocity_field(filename, output_file=None):
    """Plot velocity field"""
    print(f"Reading velocity field data from {filename}")
    # Read the data using numpy for more robust parsing
    try:
        # First, let's examine the file contents to understand its structure
        with open(filename, 'r') as f:
            first_few_lines = [next(f) for _ in range(5)]
        print(f"First few lines of {filename}:")
        for i, line in enumerate(first_few_lines):
            print(f"Line {i+1}: {line.strip()}")
        
        # Now try to load the data
        data = np.loadtxt(filename)
        print(f"Velocity field data shape: {data.shape}")
        
        # Check if the file is empty or has the wrong format
        if data.size == 0:
            print(f"Warning: {filename} appears to be empty.")
            return
            
        # Determine dimensions based on the shape of the loaded data
        ny, ncols = data.shape
        nx = ncols // 2  # Each row has pairs of u, v values
        print(f"Detected grid size: {nx}x{ny}")
        
        # Reshape the data into separate u and v arrays
        u = np.zeros((ny, nx))
        v = np.zeros((ny, nx))
        
        for j in range(ny):
            for i in range(nx):
                u[j,i] = data[j, 2*i]
                v[j,i] = data[j, 2*i+1]
        
        # Output some statistics about the velocity data
        print(f"u velocity range: {u.min()} to {u.max()}")
        print(f"v velocity range: {v.min()} to {v.max()}")
        
        # If velocities are all zero, there's an issue with the simulation output
        if np.allclose(u, 0) and np.allclose(v, 0):
            print("WARNING: All velocity values are close to zero!")
        
        # Create coordinates
        x = np.linspace(0, nx-1, nx)
        y = np.linspace(0, ny-1, ny)
        X, Y = np.meshgrid(x, y)
        
        # Calculate velocity magnitude
        speed = np.sqrt(u*u + v*v)
        print(f"Velocity magnitude range: {speed.min()} to {speed.max()}")
        
        # Plot
        plt.figure(figsize=(12, 8))
        
        # Try different visualization approaches
        try:
            # Try streamplot for well-behaved velocity fields
            plt.streamplot(X, Y, u, v, density=1.5, color=speed, cmap='viridis')
            plt.colorbar(label='Velocity magnitude')
        except Exception as stream_err:
            print(f"Streamplot failed: {stream_err}. Trying quiver plot instead.")
            # Fallback to quiver plot which is more robust
            skip = max(1, min(nx, ny) // 25)  # Skip some points to avoid crowding
            plt.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
                      u[::skip, ::skip], v[::skip, ::skip], 
                      speed[::skip, ::skip], cmap='viridis')
            plt.colorbar(label='Velocity magnitude')
        
        plt.title('Velocity Field')
        plt.xlabel('x')
        plt.ylabel('y')
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved velocity field plot to {output_file}")
        else:
            plt.show()
        plt.close()
        
    except Exception as e:
        print(f"Error plotting velocity field from {filename}: {e}")
        import traceback
        traceback.print_exc()

def create_animation_script():
    """Create a script to make an animation of the velocity profile evolution"""
    print("Creating animation script...")
    animation_script = """#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import sys

# Get all profile files
profile_files = sorted(glob.glob('output/profile_*.dat'), 
                      key=lambda x: int(x.split('_')[1].split('.')[0]))

# Print debugging info
print(f"Found {len(profile_files)} profile files")
if profile_files:
    print(f"First file: {profile_files[0]}")
    print(f"Last file: {profile_files[-1]}")

if not profile_files:
    print("No profile files found! Make sure simulation has run correctly.")
    sys.exit(1)

# Check if we can read the first file
try:
    print(f"Attempting to read {profile_files[0]}...")
    test_data = np.loadtxt(profile_files[0])
    print(f"Successfully read data with shape {test_data.shape}")
except Exception as e:
    print(f"ERROR: Failed to read first profile file: {e}")
    sys.exit(1)

fig, ax = plt.subplots(figsize=(10, 6))

def init():
    ax.set_xlabel('Velocity u$_x$')
    ax.set_ylabel('Channel height y')
    ax.set_title('Poiseuille Flow Profile Evolution')
    ax.grid(True)
    return []

def animate(i):
    ax.clear()
    print(f"Animating frame {i}/{len(profile_files)}")
    data = np.loadtxt(profile_files[i])
    y = data[:,0]
    u_sim = data[:,1]
    u_analytical = data[:,2]
    
    ax.plot(u_sim, y, 'o-', label='Simulation')
    ax.plot(u_analytical, y, '--', label='Analytical')
    ax.set_xlabel('Velocity u$_x$')
    ax.set_ylabel('Channel height y')
    ax.set_title(f'Poiseuille Flow Profile Evolution - Step {i*100}')
    ax.legend()
    ax.grid(True)
    
    return []

# Create animation with a smaller number of frames if there are too many
num_frames = min(len(profile_files), 100)
frame_step = max(1, len(profile_files) // num_frames)
frame_indices = range(0, len(profile_files), frame_step)

print(f"Creating animation with {len(frame_indices)} frames...")
ani = animation.FuncAnimation(fig, animate, frames=frame_indices,
                              init_func=init, blit=True, interval=200)

# Check for ffmpeg
import subprocess
try:
    subprocess.run(['ffmpeg', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    print("ffmpeg is available")
    # Save animation
    print("Saving animation to velocity_evolution.mp4...")
    ani.save('velocity_evolution.mp4', writer='ffmpeg', dpi=200)
    print("Animation saved successfully!")
except Exception as e:
    print(f"Could not save animation: {e}")
    print("Falling back to showing the animation instead")
    plt.show()

plt.close()
"""
    with open('animate_profiles.py', 'w') as f:
        f.write(animation_script)
    print("Animation script created as 'animate_profiles.py'")

def plot_early_timesteps():
    """Plot results from early time steps where there was actual flow before instability"""
    print("\nAnalyzing early time steps:")
    print("================================")
    
    # Create directory for early results
    os.makedirs('figures/early', exist_ok=True)
    
    # Check steps 100, 200, 300 and 400
    early_steps = [100, 200, 300]
    
    for step in early_steps:
        vel_file = f'output/vel_{step}.dat'
        profile_file = f'output/profile_{step}.dat'
        
        if os.path.exists(vel_file):
            print(f"\nAnalyzing step {step}:")
            print("-----------------------")
            # Plot velocity field
            plot_velocity_field(vel_file, f'figures/early/velocity_field_{step}.png')
            
        if os.path.exists(profile_file):
            # Plot velocity profile
            plot_velocity_profile(profile_file, f'figures/early/profile_{step}.png')
    
    print("\nEarly time step analysis complete.")
    print("Check figures/early directory for the results.")

def main():
    # Create output directory for figures
    os.makedirs('figures', exist_ok=True)
    
    # Check if output directory exists
    if not os.path.exists('output'):
        print("Error: 'output' directory not found. Make sure the simulation has run first.")
        return
    
    # Check simulation output files
    print("\nChecking simulation output files:")
    print("================================")
    for pattern in ['output/rho_*.dat', 'output/vel_*.dat', 'output/profile_*.dat']:
        files = glob.glob(pattern)
        print(f"{pattern}: Found {len(files)} files")
        if files:
            print(f"  Example: {files[0]}")
            # Check file size
            size = os.path.getsize(files[0])
            print(f"  File size: {size} bytes")
            # Check if empty
            if size == 0:
                print("  WARNING: File appears to be empty!")
    
    # Plot final velocity profile
    print("\nGenerating velocity profile:")
    print("================================")
    if os.path.exists('output/final_profile.dat'):
        plot_velocity_profile('output/final_profile.dat', 'figures/final_profile.png')
    else:
        print("Warning: output/final_profile.dat not found.")
    
    # Plot velocity field (using a known good step)
    print("\nGenerating velocity field from step 200 (before instability):")
    print("================================")
    if os.path.exists('output/vel_200.dat'):
        plot_velocity_field('output/vel_200.dat', 'figures/velocity_field_step200.png')
    else:
        # Fall back to the final step
        vel_files = sorted(glob.glob('output/vel_*.dat'), 
                            key=lambda x: int(x.split('_')[1].split('.')[0]))
        if vel_files:
            last_vel_file = vel_files[-1]
            plot_velocity_field(last_vel_file, 'figures/velocity_field.png')
        else:
            print("Warning: No velocity files found in output directory.")
    
    # Create animation script
    print("\nCreating animation script:")
    print("================================")
    create_animation_script()
    
    # Analyze early timesteps
    plot_early_timesteps()
    
    print("\nSummary:")
    print("================================")
    print("Plots generated in the 'figures' directory.")
    print("Early-time-step plots are in 'figures/early' directory.")
    print("To create an animation, run: python animate_profiles.py")
    print("If you encounter issues, check if ffmpeg is installed on your system.")

if __name__ == "__main__":
    print("Starting visualization...")
    main() 