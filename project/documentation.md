# Entropic Lattice Boltzmann Method Simulation

This documentation provides a comprehensive overview of the Entropic Lattice Boltzmann Method (ELBM) simulation codebase. The simulation is designed to model fluid flow using the lattice Boltzmann method with entropic stabilization, which ensures numerical stability by enforcing the H-theorem.

## Table of Contents

1. [Project Structure](#project-structure)
2. [Core Components](#core-components)
3. [Collision Models](#collision-models)
4. [Test Cases](#test-cases)
5. [Boundary Conditions](#boundary-conditions)
6. [Visualization](#visualization)
7. [Building and Running](#building-and-running)
8. [Testing](#testing)

## Project Structure

The project is organized into the following directories and files:

```
project/
├── include/                 # Header files
│   ├── boundary_conditions.h
│   ├── common_types.h
│   ├── data_writer.h
│   ├── entropic_collision.h
│   ├── entropic_equilibrium.h
│   ├── lattice.h
│   ├── node_data.h
│   ├── numerical_solvers.h
│   ├── simulator.h
│   ├── test_cases.h
│   ├── trt_collision.h
│   └── vtk_writer.h
├── src/                     # Implementation files
│   ├── boundary_conditions.cpp
│   ├── data_writer.cpp
│   ├── entropic_collision.cpp
│   ├── entropic_equilibrium.cpp
│   ├── lattice.cpp
│   ├── main.cpp
│   ├── node_data.cpp
│   ├── numerical_solvers.cpp
│   ├── simulator.cpp
│   ├── simulator_geometries.cpp
│   ├── test_cases.cpp
│   ├── trt_collision.cpp
│   └── vtk_writer.cpp
├── figures/                 # Output directory for visualization
├── results/                 # Output directory for simulation results
├── Makefile                 # Build system
├── report.tex               # LaTeX report
├── visualize_results.py     # Python visualization script
└── documentation.md         # This documentation file
```

## Core Components

### Common Types (`common_types.h`)

Defines common types used throughout the codebase:

- `Real`: Floating-point type (double precision)
- `Vector2D`: 2D vector structure with x and y components
- `REAL_EPSILON`: Small positive value to prevent division by zero

### Lattice (`lattice.h`, `lattice.cpp`)

Defines the D2Q9 lattice structure used in the simulation:

- 9 discrete velocities in 2D
- Weights for each direction
- Speed of sound
- Opposite directions for each velocity

### Node Data (`node_data.h`, `node_data.cpp`)

Represents the state of a single grid node:

- Distribution functions (`f`)
- Equilibrium distribution functions (`f_eq`)
- Post-streaming distribution functions (`f_new`)
- Macroscopic variables: density (`rho`) and velocity (`u`)
- Fluid/solid flag (`is_fluid`)
- Methods for computing macroscopic variables and equilibrium distributions

### Simulator (`simulator.h`, `simulator.cpp`, `simulator_geometries.cpp`)

Central class that orchestrates the simulation:

- Initializes the grid
- Sets up boundary conditions
- Performs time stepping (collision, streaming, boundary conditions)
- Tracks convergence
- Writes output data

### Data Writer (`data_writer.h`, `data_writer.cpp`)

Handles output of simulation results:

- Writes grid data to .dat files
- Formats data for visualization

## Collision Models

The simulation supports two collision models:

### Entropic Collision (`entropic_collision.h`, `entropic_collision.cpp`)

Implements the entropic collision operator:

- Ensures numerical stability by enforcing the H-theorem
- Uses a variable relaxation parameter α determined by solving a nonlinear equation
- Requires more computational resources but provides better stability

### Two-Relaxation-Time (TRT) Collision (`trt_collision.h`, `trt_collision.cpp`)

Implements the TRT collision operator:

- Separates relaxation of symmetric (even) and antisymmetric (odd) parts
- Uses two relaxation times: τ+ (related to viscosity) and τ- (set using magic parameter)
- Provides better accuracy for boundary conditions
- Computationally more efficient than entropic collision

## Test Cases

The simulation includes several test cases to benchmark the collision models:

### Elliptical Flow (`test_cases.cpp`)

- Flow around an elliptical obstacle
- Classic benchmark for external flows
- Demonstrates vortex shedding at moderate Reynolds numbers

### Lid-Driven Cavity (`test_cases.cpp`)

- Flow in a square cavity with a moving top wall
- Standard benchmark for internal flows
- Tests accuracy of boundary conditions

### Channel Flow with Obstacle (`test_cases.cpp`)

- Flow in a channel with a circular obstacle
- Similar to elliptical flow but with different geometry
- Tests the ability to handle curved boundaries

### Backward-Facing Step (`test_cases.cpp`)

- Flow over a backward-facing step
- Tests flow separation and reattachment
- Demonstrates recirculation zones

### Taylor-Green Vortex (`test_cases.cpp`)

- Decaying vortex flow with analytical solution
- Tests accuracy of the numerical scheme
- Uses periodic boundary conditions

### Poiseuille Flow (`test_cases.cpp`)

- Pressure-driven flow in a channel
- Has analytical solution for comparison
- Tests pressure boundary conditions

## Boundary Conditions

The simulation supports various boundary conditions:

### Bounce-Back (`boundary_conditions.cpp`)

- No-slip boundary condition for walls
- Simple and robust implementation
- Second-order accurate for straight walls

### Zou/He Velocity Boundary (`boundary_conditions.cpp`)

- Implements fixed velocity boundary condition
- Used for inlet and lid-driven cavity
- Maintains mass conservation

### Pressure Boundary (`boundary_conditions.cpp`)

- Implements fixed pressure (density) boundary condition
- Used for Poiseuille flow
- Allows pressure-driven flows

### Periodic Boundary (`boundary_conditions.cpp`)

- Connects opposite sides of the domain
- Used for Taylor-Green vortex
- Simulates infinite domain

## Visualization

The simulation includes a Python script for visualizing results:

### `visualize_results.py`

- Reads .dat files produced by the simulation
- Creates plots of velocity magnitude, vorticity, and streamlines
- Generates animations of the flow evolution
- Plots convergence history

## Building and Running

### Building the Code

The code can be built using the provided Makefile:

```bash
make clean   # Clean previous build
make         # Build the code
```

### Running the Simulation

The simulation can be run with various command-line options:

```bash
./entropic_lbm                      # Run with default parameters
./entropic_lbm --test ellipse_flow  # Run elliptical flow test case
./entropic_lbm --test lid_driven_cavity  # Run lid-driven cavity test case
./entropic_lbm --nx 200 --ny 100    # Set grid size
./entropic_lbm --re 100             # Set Reynolds number
./entropic_lbm --trt                # Use TRT collision model
./entropic_lbm --magic 0.25         # Set magic parameter for TRT
./entropic_lbm --list-tests         # List available test cases
./entropic_lbm --help               # Show help message
```

### Visualizing Results

The visualization script can be run with various options:

```bash
python visualize_results.py                 # Visualize with default options
python visualize_results.py --type velocity  # Visualize velocity field
python visualize_results.py --type vorticity # Visualize vorticity field
python visualize_results.py --animate        # Create animations
python visualize_results.py --convergence    # Plot convergence history
python visualize_results.py --step 1000      # Visualize specific time step
python visualize_results.py --test ellipse_flow_entropic  # Specify test case for file naming
```

The `--test` option allows you to specify a test case name for file naming, which is useful when comparing different test cases or collision models.

## Testing

The code includes a comprehensive test suite to ensure correctness:

### Test Cases

Each test case can be run individually using the provided Makefile targets:

```bash
make test-ellipse      # Run elliptical flow test case with both collision models
make test-cavity       # Run lid-driven cavity test case with both collision models
make test-channel      # Run channel flow with obstacle test case with both collision models
make test-step         # Run backward-facing step test case with both collision models
make test-vortex       # Run Taylor-Green vortex test case with both collision models
make test-poiseuille   # Run Poiseuille flow test case with both collision models
```

Alternatively, you can run individual test cases directly:

```bash
./entropic_lbm --test ellipse_flow
./entropic_lbm --test lid_driven_cavity
./entropic_lbm --test channel_obstacle
./entropic_lbm --test backward_facing_step
./entropic_lbm --test taylor_green_vortex
./entropic_lbm --test poiseuille_flow
```

### Collision Model Comparison

The entropic and TRT collision models can be compared by running the same test case with different models:

```bash
./entropic_lbm --test ellipse_flow           # Use entropic collision
./entropic_lbm --test ellipse_flow --trt     # Use TRT collision
```

You can also use the dedicated comparison target:

```bash
make compare-models    # Run elliptical flow with both models and save results in separate directories
```

### Convergence Testing

The convergence of the simulation can be monitored by plotting the convergence history:

```bash
python visualize_results.py --convergence
```

This will create a plot of the maximum velocity change over time, which should decrease as the simulation converges.

### Cleaning Up

To clean up the build and results:

```bash
make clean             # Remove object files, executable, and results
make clean-results     # Remove only the results, keeping the executable
```
