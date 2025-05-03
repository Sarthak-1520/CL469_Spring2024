#!/bin/bash

# Script to run LBM simulation for different parameter sets
# required for Assignment Two report figures.

# Exit immediately if a command exits with a non-zero status.
set -e

# Base directory where the script is located (assignment_two)
BASE_DIR=$(pwd)
PARAM_FILE="${BASE_DIR}/parameters.dat"
DATA_BASE_DIR="${BASE_DIR}/data"
EXE="${BASE_DIR}/build/lbm_poiseuille"

# Function to modify parameters.dat using sed
# Usage: modify_param parameter_name new_value
modify_param() {
    local param_key=$1
    local new_val=$2
    # Use sed to find the line starting with the key (allowing spaces) and replace the value after '='
    # Use a temporary file for macOS compatibility with sed -i
    sed -i.bak "s/^[[:space:]]*${param_key}[[:space:]]*=.*$/${param_key} = ${new_val}/" "${PARAM_FILE}"
    rm -f "${PARAM_FILE}.bak" # Clean up backup file
    echo "Set ${param_key} = ${new_val} in ${PARAM_FILE}"
}

# Function to run simulation and store results
# Usage: run_and_store H tau collision_op [tau_minus]
run_and_store() {
    local H_val=$1
    local tau_val=$2
    local coll_op=$3
    local tau_m_val=${4:-1.0} # Default tau_minus if not provided

    # --- Modify parameters.dat ---
    modify_param "H" ${H_val}
    modify_param "tau" ${tau_val}
    modify_param "collision_operator" ${coll_op}
    modify_param "tau_minus" ${tau_m_val}

    # Recalculate nu and g based on H and tau (to keep Re constant for resolution study?)
    # Assignment text implies nu=0.1, g=9.203e-6 for H=11, tau=0.8 base case.
    # For resolution study (Sec 5), H changes. Re = um*H/nu. Need to keep Re same.
    # um = g*H^2/(8*nu). Re = g*H^3/(8*nu^2). Let Re_target = 0.1531.
    # Let's fix nu = (tau-0.5)/3. g = sqrt(Re_target * 8 * nu^2 / H^3)? No, assignment says vary H, not g/nu explicitly for force error plot.
    # Let's assume for Sec 5 (Force Error vs H), we just change H and keep other base parameters (tau=0.8, nu=0.1, g=9.203e-6) fixed.
    # This means Re WILL change, but assignment asks to plot error vs H anyway.

    # For slip velocity study (Sec 6), tau changes. Keep H=11, g=9.203e-6 fixed.
    # nu = (tau-0.5)/3. Need to recalculate nu based on tau.
    if [[ "${coll_op}" == "BGK" || "${coll_op}" == "TRT" ]]; then
        local nu_val=$(echo "(${tau_val} - 0.5) / 3.0" | bc -l)
        modify_param "nu" ${nu_val}
        echo "Calculated nu = ${nu_val} for tau = ${tau_val}"
    fi

    # --- Define output directory ---
    local output_subdir=""
    # Naming convention based on plot_report_figures.py expectations
    if [[ "${coll_op}" == "RES" ]]; then # Special marker for resolution study
        output_subdir="data_H${H_val}"
        coll_op="BGK" # Use BGK for resolution study unless specified otherwise
        modify_param "collision_operator" ${coll_op}
    elif [[ "${coll_op}" == "BGK" ]]; then
         output_subdir="data_BGK_tau${tau_val}"
    elif [[ "${coll_op}" == "TRT" ]]; then
         output_subdir="data_TRT_tau${tau_val}" # Using tau+ (which is tau) for naming
    else
        echo "Error: Unknown collision operator type for directory naming: ${coll_op}"
        exit 1
    fi
    local output_dir="${DATA_BASE_DIR}/${output_subdir}"

    echo "--------------------------------------------------"
    echo "Running for: H=${H_val}, tau=${tau_val}, Operator=${coll_op}, nu=${nu_val:-$(grep '^nu' ${PARAM_FILE} | awk '{print $3}')}"
    echo "Output Directory: ${output_dir}"
    echo "--------------------------------------------------"

    # --- Clean old data, Run Simulation ---
    rm -rf "${DATA_BASE_DIR}"/*.dat "${DATA_BASE_DIR}"/results_summary.txt # Clear base data dir
    make run # Runs the simulation which outputs to DATA_BASE_DIR

    # --- Store results ---
    mkdir -p "${output_dir}"
    echo "Moving results to ${output_dir}"
    mv "${DATA_BASE_DIR}"/*.dat "${DATA_BASE_DIR}"/results_summary.txt "${output_dir}/"

    # --- Restore base parameters? (Optional, but safer) ---
    # modify_param H 11
    # modify_param tau 0.8
    # modify_param collision_operator BGK
    # modify_param nu 0.1
    # modify_param g 9.203e-6
}

# === Parameter Sets ===

# --- Resolution Study (Report Section 5) ---
# Vary H, keep tau=0.8, nu=0.1, g=9.203e-6, use BGK
H_values=(5 10 20 40)
ref_tau=0.8
ref_coll_op="RES" # Use special marker

echo "=== Starting Resolution Study ==="
for H in "${H_values[@]}"; do
    run_and_store ${H} ${ref_tau} ${ref_coll_op}
done

# --- Slip Velocity Study (Report Section 6) ---
# Vary tau, keep H=11, g=9.203e-6. Compare BGK vs TRT.
# For TRT, use tau_minus = 1.0 (default) or maybe the magic parameter lambda=1/4 value?
# lambda = (1/tau_m - 0.5)*(1/tau_p - 0.5) = 1/4.
# tau_p = nu/cs2 + 0.5 = (tau-0.5)/cs2/3 + 0.5 = tau.
# (1/tau_m - 0.5)*(1/tau - 0.5) = 1/4.
# Let's just use tau_minus=1.0 for simplicity as done in Parameters.cpp default.
tau_slip_values=(0.8 2.0 5.0)
ref_H=11
tau_m_slip=1.0

echo "=== Starting Slip Velocity Study (BGK) ==="
for tau in "${tau_slip_values[@]}"; do
    run_and_store ${ref_H} ${tau} "BGK"
done

echo "=== Starting Slip Velocity Study (TRT) ==="
for tau in "${tau_slip_values[@]}"; do
    run_and_store ${ref_H} ${tau} "TRT" ${tau_m_slip}
done

# --- Restore parameters.dat to original state (optional) ---
echo "Restoring parameters.dat to default H=11, tau=0.8, BGK..."
modify_param "H" 11
modify_param "tau" 0.8
modify_param "collision_operator" BGK
modify_param "nu" 0.1
modify_param "g" 9.203e-6
modify_param "tau_minus" 1.0

echo "=== Parameter study complete ==="

# --- Copy Reference Case Data Back for Basic Plots ---
echo "=== Copying reference case data (BGK, tau=0.8) for basic plots ==="
mkdir -p "${DATA_BASE_DIR}"
cp "${DATA_BASE_DIR}/data_BGK_tau0.8/"* "${DATA_BASE_DIR}/"

# --- Generate Plots ---
echo "=== Generating final visualizations ==="
make visualize

echo "=== Script finished ===" 