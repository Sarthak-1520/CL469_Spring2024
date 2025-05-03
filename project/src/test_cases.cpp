#include "test_cases.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

// Factory function to create test cases
std::unique_ptr<TestCase> createTestCase(const std::string& name) {
    if (name == "ellipse_flow") {
        return std::make_unique<EllipseFlowTestCase>();
    } else if (name == "lid_driven_cavity") {
        return std::make_unique<LidDrivenCavityTestCase>();
    } else if (name == "channel_obstacle") {
        return std::make_unique<ChannelObstacleTestCase>();
    } else if (name == "backward_facing_step") {
        return std::make_unique<BackwardFacingStepTestCase>();
    } else if (name == "taylor_green_vortex") {
        return std::make_unique<TaylorGreenVortexTestCase>();
    } else if (name == "poiseuille_flow") {
        return std::make_unique<PoiseuilleFlowTestCase>();
    } else {
        throw std::invalid_argument("Unknown test case: " + name);
    }
}

// Get a list of all available test cases
std::vector<std::string> getAvailableTestCases() {
    return {
        "ellipse_flow",
        "lid_driven_cavity",
        "channel_obstacle",
        "backward_facing_step",
        "taylor_green_vortex",
        "poiseuille_flow"
    };
}

// ===== Ellipse Flow Test Case =====
std::unique_ptr<Simulator> EllipseFlowTestCase::setup(
    int nx, int ny, Real Re, bool use_trt, Real magic_param) {

    // Ellipse parameters
    Real ellipse_cx = nx / 4.0;      // Center X (place upstream)
    Real ellipse_cy = ny / 2.0;      // Center Y
    Real ellipse_a = ny / 8.0;       // Semi-axis X
    Real ellipse_b = ny / 8.0;       // Semi-axis Y

    // Inlet velocity based on low Mach number
    Real cs = 1.0 / std::sqrt(3.0);
    Real target_Ma = 0.05;
    Real U_in = target_Ma * cs; // Inlet velocity

    // Characteristic length L = 2 * ellipse_b (minor diameter)
    Real characteristic_length = 2.0 * ellipse_b;
    // viscosity = U * L / Re
    Real viscosity = U_in * characteristic_length / Re;

    int total_steps = 1000;    // Total simulation steps (reduced for testing)
    int output_freq = 200;     // Frequency of writing output files

    std::cout << "--- Ellipse Flow Parameters ---" << std::endl;
    std::cout << "Center: (" << ellipse_cx << ", " << ellipse_cy << ")" << std::endl;
    std::cout << "Semi-axes: a=" << ellipse_a << ", b=" << ellipse_b << std::endl;
    std::cout << "Inlet U: " << U_in << " (Ma=" << U_in/cs << ")" << std::endl;
    std::cout << "Re = " << Re << ", Viscosity = " << viscosity << std::endl;
    std::cout << "-------------------------" << std::endl;

    // Create simulator with selected collision model
    auto sim = std::make_unique<Simulator>(
        nx, ny, viscosity, U_in, total_steps, output_freq,
        use_trt ? Simulator::CollisionModel::TRT : Simulator::CollisionModel::ENTROPIC,
        magic_param,
        GeometryType::ELLIPSE,
        ellipse_cx, ellipse_cy, ellipse_a, ellipse_b
    );

    return sim;
}

// ===== Lid-Driven Cavity Test Case =====
std::unique_ptr<Simulator> LidDrivenCavityTestCase::setup(
    int nx, int ny, Real Re, bool use_trt, Real magic_param) {

    // Lid velocity based on low Mach number
    Real cs = 1.0 / std::sqrt(3.0);
    Real target_Ma = 0.1;
    Real U_lid = target_Ma * cs; // Lid velocity

    // Characteristic length L = cavity height
    Real characteristic_length = static_cast<Real>(ny - 1);
    // viscosity = U * L / Re
    Real viscosity = U_lid * characteristic_length / Re;

    int total_steps = 1000;    // Total simulation steps (reduced for testing)
    int output_freq = 200;     // Frequency of writing output files

    std::cout << "--- Lid-Driven Cavity Parameters ---" << std::endl;
    std::cout << "Cavity size: " << nx << " x " << ny << std::endl;
    std::cout << "Lid velocity: " << U_lid << " (Ma=" << U_lid/cs << ")" << std::endl;
    std::cout << "Re = " << Re << ", Viscosity = " << viscosity << std::endl;
    std::cout << "-------------------------" << std::endl;

    // Create simulator with selected collision model
    auto sim = std::make_unique<Simulator>(
        nx, ny, viscosity, U_lid, total_steps, output_freq,
        use_trt ? Simulator::CollisionModel::TRT : Simulator::CollisionModel::ENTROPIC,
        magic_param,
        GeometryType::LID_DRIVEN_CAVITY
    );

    return sim;
}

// ===== Channel Flow with Obstacle Test Case =====
std::unique_ptr<Simulator> ChannelObstacleTestCase::setup(
    int nx, int ny, Real Re, bool use_trt, Real magic_param) {

    // Obstacle parameters
    Real obstacle_cx = nx / 3.0;     // Center X
    Real obstacle_cy = ny / 2.0;     // Center Y
    Real obstacle_r = ny / 8.0;      // Radius

    // Inlet velocity based on low Mach number
    Real cs = 1.0 / std::sqrt(3.0);
    Real target_Ma = 0.05;
    Real U_in = target_Ma * cs; // Inlet velocity

    // Characteristic length L = 2 * obstacle_r (diameter)
    Real characteristic_length = 2.0 * obstacle_r;
    // viscosity = U * L / Re
    Real viscosity = U_in * characteristic_length / Re;

    int total_steps = 1000;    // Total simulation steps (reduced for testing)
    int output_freq = 200;     // Frequency of writing output files

    std::cout << "--- Channel Flow with Obstacle Parameters ---" << std::endl;
    std::cout << "Obstacle center: (" << obstacle_cx << ", " << obstacle_cy << ")" << std::endl;
    std::cout << "Obstacle radius: " << obstacle_r << std::endl;
    std::cout << "Inlet U: " << U_in << " (Ma=" << U_in/cs << ")" << std::endl;
    std::cout << "Re = " << Re << ", Viscosity = " << viscosity << std::endl;
    std::cout << "-------------------------" << std::endl;

    // Create simulator with selected collision model
    auto sim = std::make_unique<Simulator>(
        nx, ny, viscosity, U_in, total_steps, output_freq,
        use_trt ? Simulator::CollisionModel::TRT : Simulator::CollisionModel::ENTROPIC,
        magic_param,
        GeometryType::CHANNEL_OBSTACLE,
        obstacle_cx, obstacle_cy, obstacle_r, 0.0 // Last parameter unused for circle
    );

    return sim;
}

// ===== Backward-Facing Step Test Case =====
std::unique_ptr<Simulator> BackwardFacingStepTestCase::setup(
    int nx, int ny, Real Re, bool use_trt, Real magic_param) {

    // Step parameters
    int step_height = ny / 2;
    int step_length = nx / 6;

    // Inlet velocity based on low Mach number
    Real cs = 1.0 / std::sqrt(3.0);
    Real target_Ma = 0.05;
    Real U_in = target_Ma * cs; // Inlet velocity

    // Characteristic length L = channel height after step
    Real characteristic_length = static_cast<Real>(ny);
    // viscosity = U * L / Re
    Real viscosity = U_in * characteristic_length / Re;

    int total_steps = 1000;    // Total simulation steps (reduced for testing)
    int output_freq = 200;     // Frequency of writing output files

    std::cout << "--- Backward-Facing Step Parameters ---" << std::endl;
    std::cout << "Step height: " << step_height << std::endl;
    std::cout << "Step length: " << step_length << std::endl;
    std::cout << "Inlet U: " << U_in << " (Ma=" << U_in/cs << ")" << std::endl;
    std::cout << "Re = " << Re << ", Viscosity = " << viscosity << std::endl;
    std::cout << "-------------------------" << std::endl;

    // Create simulator with selected collision model
    auto sim = std::make_unique<Simulator>(
        nx, ny, viscosity, U_in, total_steps, output_freq,
        use_trt ? Simulator::CollisionModel::TRT : Simulator::CollisionModel::ENTROPIC,
        magic_param,
        GeometryType::BACKWARD_FACING_STEP,
        step_length, step_height, 0.0, 0.0 // Last two parameters unused
    );

    return sim;
}

// ===== Taylor-Green Vortex Test Case =====
std::unique_ptr<Simulator> TaylorGreenVortexTestCase::setup(
    int nx, int ny, Real Re, bool use_trt, Real magic_param) {

    // Taylor-Green vortex parameters
    Real U0 = 0.05; // Maximum velocity

    // Characteristic length L = domain size
    Real characteristic_length = static_cast<Real>(nx);
    // viscosity = U * L / Re
    Real viscosity = U0 * characteristic_length / Re;

    int total_steps = 1000;    // Total simulation steps (reduced for testing)
    int output_freq = 200;     // Frequency of writing output files

    std::cout << "--- Taylor-Green Vortex Parameters ---" << std::endl;
    std::cout << "Domain size: " << nx << " x " << ny << std::endl;
    std::cout << "Maximum velocity: " << U0 << std::endl;
    std::cout << "Re = " << Re << ", Viscosity = " << viscosity << std::endl;
    std::cout << "-------------------------" << std::endl;

    // Create simulator with selected collision model
    auto sim = std::make_unique<Simulator>(
        nx, ny, viscosity, U0, total_steps, output_freq,
        use_trt ? Simulator::CollisionModel::TRT : Simulator::CollisionModel::ENTROPIC,
        magic_param,
        GeometryType::TAYLOR_GREEN_VORTEX
    );

    return sim;
}

// ===== Poiseuille Flow Test Case =====
std::unique_ptr<Simulator> PoiseuilleFlowTestCase::setup(
    int nx, int ny, Real Re, bool use_trt, Real magic_param) {

    // Poiseuille flow parameters
    Real U_max = 0.05; // Maximum velocity at center

    // Characteristic length L = channel height
    Real characteristic_length = static_cast<Real>(ny - 1);
    // viscosity = U * L / Re
    Real viscosity = U_max * characteristic_length / Re;

    int total_steps = 1000;    // Total simulation steps (reduced for testing)
    int output_freq = 200;     // Frequency of writing output files

    std::cout << "--- Poiseuille Flow Parameters ---" << std::endl;
    std::cout << "Channel size: " << nx << " x " << ny << std::endl;
    std::cout << "Maximum velocity: " << U_max << std::endl;
    std::cout << "Re = " << Re << ", Viscosity = " << viscosity << std::endl;
    std::cout << "-------------------------" << std::endl;

    // Create simulator with selected collision model
    auto sim = std::make_unique<Simulator>(
        nx, ny, viscosity, U_max, total_steps, output_freq,
        use_trt ? Simulator::CollisionModel::TRT : Simulator::CollisionModel::ENTROPIC,
        magic_param,
        GeometryType::POISEUILLE_FLOW
    );

    return sim;
}
