# Entropic Lattice Boltzmann Method Simulation

This project implements an Entropic Lattice Boltzmann Method (ELBM) for simulating various fluid flow problems. The entropic approach ensures numerical stability by enforcing the H-theorem through a variable relaxation parameter, making it particularly suitable for high Reynolds number flows and complex geometries.

## Building the Project

To build the project, run:

```bash
make
```

This will compile the source code and create the executable `entropic_lbm`.

## Running Simulations

### Default Simulation

To run the default simulation (elliptical flow) and generate visualizations:

```bash
make run
```

### Data Generation Only

To run the default simulation and generate data files only (without visualization):

```bash
make run-data-only
```

### Test Cases

The project includes several benchmark test cases:

1. **Elliptical Flow**:
   - With visualization: `make test-ellipse`
   - Data only: `make data-ellipse`

2. **Lid-Driven Cavity Flow**:
   - With visualization: `make test-cavity`
   - Data only: `make data-cavity`

3. **Channel Flow with Circular Obstacle**:
   - With visualization: `make test-channel`
   - Data only: `make data-channel`

4. **Backward-Facing Step Flow**:
   - With visualization: `make test-step`
   - Data only: `make data-step`

5. **Taylor-Green Vortex**:
   - With visualization: `make test-vortex`
   - Data only: `make data-vortex`

6. **Poiseuille Flow**:
   - With visualization: `make test-poiseuille`
   - Data only: `make data-poiseuille`

### Generate All Data Files

To generate data files for all test cases:

```bash
make data-all
```

## Command Line Options

The `entropic_lbm` executable supports several command line options:

```bash
./entropic_lbm [options]
```

Options:
- `--test TESTNAME`: Run a specific test case (default: ellipse_flow)
- `--nx VALUE`: Grid size in x direction (default: 400)
- `--ny VALUE`: Grid size in y direction (default: 100)
- `--re VALUE`: Reynolds number (default: 100)
- `--trt`: Use TRT collision model (default: entropic)
- `--magic VALUE`: Magic parameter for TRT (default: 0.25)
- `--list-tests`: List available test cases
- `--help`: Show help message

## Visualization

The Python script `visualize_results.py` is used to create visualizations from the simulation data:

```bash
python visualize_results.py [options]
```

Options:
- `--type {velocity,vorticity,streamlines,all}`: Type of plot to generate
- `--animate`: Create animation
- `--convergence`: Plot convergence history
- `--step STEP`: Specific time step to visualize
- `--test TEST`: Test case name for file naming

## Cleaning Up

To clean up object files and the executable:

```bash
make clean
```

To clean up results only:

```bash
make clean-results
```

## Output Files

The simulation generates several types of output files:

1. **Data Files**: Located in the `results` directory with the pattern `lbm_output_*.dat`
2. **Convergence Data**: Located at `results/convergence.dat`
3. **Visualization Images**: Located in the `figures` directory

## Collision Models

The simulation supports two collision models:

1. **Entropic Collision Model**: Ensures numerical stability by enforcing the H-theorem
2. **Two-Relaxation-Time (TRT) Model**: Offers a compromise between simplicity and stability

## Project Structure

```
project/
├── include/                 # Header files
├── src/                     # Implementation files
├── results/                 # Output data files
├── figures/                 # Visualization images
├── Makefile                 # Build system
├── visualize_results.py     # Visualization script
└── README.md                # This file
```