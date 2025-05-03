# CL469 Assignment Two: LBM Simulation of Poiseuille Flow

This project simulates 2D gravity-driven Poiseuille flow between parallel plates using the Lattice Boltzmann Method (LBM) with the D2Q9 model, following the requirements of Assignment Two for CL469.

## Problem Description

The simulation models water flow in a channel of physical length Lp = 0.1 m and width Hp = 5e-5 m under gravity gx = 9.8 m/s^2. The kinematic viscosity is nu_p = 1e-6 m^2/s and density rho_p = 1000 kg/m^3.

The simulation uses lattice units, matching the physical Reynolds number Re = 0.1531 (for the base case). Periodic boundary conditions are applied in the x-direction, and no-slip (bounce-back) conditions are applied at the top and bottom walls. Both BGK and TRT collision operators are implemented, along with Guo forcing for the body force.

## Project Structure

```
assignment_two/
├── Makefile                  # Build script
├── README.md                 # This file
├── parameters.dat            # Simulation parameters (edited by study script)
├── report.tex                # LaTeX report source
├── report.pdf                # Compiled report (generated)
├── run_parameter_study.sh    # Script to automate simulation runs
├── include/                  # Header files (.h)
│   ├── Lattice.h
│   ├── Parameters.h
│   └── Simulation.h
├── src/                      # Source files (.cpp)
│   ├── Parameters.cpp
│   ├── Simulation.cpp
│   └── main.cpp
├── build/                    # Compiled object files and executable (generated)
├── data/                     # Simulation output data (generated)
│   ├── data_H5/              # Data for H=5 resolution study
│   ├── data_H10/
│   ├── data_H20/
│   ├── data_H40/
│   ├── data_BGK_tau0.8/      # Data for BGK slip study (tau=0.8)
│   ├── data_BGK_tau2.0/
│   ├── data_BGK_tau5.0/
│   ├── data_TRT_tau0.8/      # Data for TRT slip study (tau=0.8)
│   ├── data_TRT_tau2.0/
│   └── data_TRT_tau5.0/
├── figures/                  # Plots for the report (generated)
│   ├── convergence.png
│   ├── density_profile.png
│   ├── dissipation_error_vs_resolution.png
│   ├── dissipation_vs_resolution.png
│   ├── force_error_vs_resolution.png
│   ├── force_methods_comparison.png
│   ├── forces_nondim_vs_time.png
│   ├── forces_vs_time.png
│   ├── forces_y_vs_time.png
│   ├── slip_velocity_comparison.png
│   ├── uy_profile.png
│   └── velocity_profile.png
└── visualizations/           # Python scripts for plotting
    ├── plot_convergence.py
    ├── plot_forces.py
    ├── plot_report_figures.py
    └── plot_velocity.py
```

## Dependencies

- C++ compiler (g++ recommended, supporting C++11)
- Make
- Python 3 with libraries (numpy, matplotlib, pandas)
- LaTeX distribution (e.g., TeX Live) for report compilation
- `bc` command-line utility (used by `run_parameter_study.sh`)

## Build and Run Instructions

1.  **Compile the code:**
    ```bash
    make
    ```
    This creates the executable `build/lbm_poiseuille`.

2.  **Run a Single Simulation (Optional):**
    To run with parameters currently in `parameters.dat`:
    ```bash
    make run
    ```
    Output data files appear in `data/`. Basic plots (for this single run) can be generated:
    ```bash
    make visualize_single # (Requires adding this target to Makefile if desired)
    ```

3.  **Run Parameter Study and Generate All Figures:**
    This is the primary way to generate results for the report.
    *   Make the script executable (only needed once):
        ```bash
        chmod +x run_parameter_study.sh
        ```
    *   Execute the script:
        ```bash
        ./run_parameter_study.sh
        ```
    This script will:
    *   Modify `parameters.dat` for each required case (varying H, tau, operator).
    *   Run the simulation (`make run`).
    *   Move the output files from `data/` to a specific subdirectory (e.g., `data/data_H5/`, `data/data_BGK_tau2.0/`).
    *   Finally, run `make visualize`, which executes all Python plotting scripts.
        *   `plot_report_figures.py` uses the data in the subdirectories to create the comparison plots (`force_error_vs_resolution.png`, `slip_velocity_comparison.png`, `dissipation_vs_resolution.png`).
        *   The other plotting scripts (`plot_velocity.py`, etc.) will likely fail at this stage because the base `data/` directory is now empty, but this is expected.

4.  **Compile the report:**
    Make sure the figures in `figures/` look correct.
    ```bash
    make report
    ```
    This generates `report.pdf` using the figures.

## Cleaning Up

To remove generated files (build artifacts, all data, figures, report files):
```bash
make clean
```

## Report Requirements Addressed

The simulation code and visualization scripts address the points outlined in the "Report preparation" section of the assignment:

1.  **Velocity Profile:** Compared with analytical solution (Fig: `velocity_profile.png` by `plot_velocity.py`).
2.  **Wall Forces (ME):** Non-dimensional forces compared to theory (Fig: `forces_nondim_vs_time.png` by `plot_forces.py`).
3.  **Mach Number Independence:** Requires manual runs by altering grid resolution and `g` while keeping `Re` and `um` constant (not automated by `run_parameter_study.sh`).
4.  **Viscous Stress Tensor:** Verification requires adding specific output/analysis code to `Simulation.cpp` (not implemented here).
5.  **Force Method Comparison:** Relative error vs. H plotted for ME, SI, FD methods (Fig: `force_error_vs_resolution.png` by `plot_report_figures.py`).
6.  **Slip Velocity (BGK vs TRT):** Comparison for high tau values (Fig: `slip_velocity_comparison.png` by `plot_report_figures.py`).
7.  **Viscous Dissipation:** Comparison with analytical value vs. H (Fig: `dissipation_vs_resolution.png` by `plot_report_figures.py`). 