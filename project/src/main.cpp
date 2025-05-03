#include "simulator.h"
#include "test_cases.h"
#include <iostream>
#include <cmath>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <memory>

void print_usage() {
    std::cout << "Usage: entropic_lbm [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --test TESTNAME    Run a specific test case (default: ellipse_flow)" << std::endl;
    std::cout << "  --nx VALUE         Grid size in x direction (default: 400)" << std::endl;
    std::cout << "  --ny VALUE         Grid size in y direction (default: 100)" << std::endl;
    std::cout << "  --re VALUE         Reynolds number (default: 100)" << std::endl;
    std::cout << "  --trt              Use TRT collision model (default: entropic)" << std::endl;
    std::cout << "  --magic VALUE      Magic parameter for TRT (default: 0.25)" << std::endl;
    std::cout << "  --list-tests       List available test cases" << std::endl;
    std::cout << "  --help             Show this help message" << std::endl;
    std::cout << std::endl;
    std::cout << "Available test cases:" << std::endl;

    auto test_cases = getAvailableTestCases();
    for (const auto& test : test_cases) {
        std::unique_ptr<TestCase> tc = createTestCase(test);
        std::cout << "  " << test << " - " << tc->getDescription() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "Entropic Lattice Boltzmann Method Simulation" << std::endl;

    // --- Create Results Directory ---
    std::filesystem::create_directory("results");

    // --- Default Parameters ---
    int nx = 400;
    int ny = 100;
    Real Re = 100.0; // Reynolds number
    bool use_trt = false; // Use entropic model by default
    Real magic_param = 0.25; // Magic parameter for TRT
    std::string test_case_name = "ellipse_flow"; // Default test case

    // --- Parse Command Line Arguments ---
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            print_usage();
            return 0;
        } else if (arg == "--list-tests") {
            std::cout << "Available test cases:" << std::endl;
            auto test_cases = getAvailableTestCases();
            for (const auto& test : test_cases) {
                std::unique_ptr<TestCase> tc = createTestCase(test);
                std::cout << "  " << test << " - " << tc->getDescription() << std::endl;
            }
            return 0;
        } else if (arg == "--test" && i + 1 < argc) {
            test_case_name = argv[++i];
        } else if (arg == "--nx" && i + 1 < argc) {
            nx = std::stoi(argv[++i]);
        } else if (arg == "--ny" && i + 1 < argc) {
            ny = std::stoi(argv[++i]);
        } else if (arg == "--re" && i + 1 < argc) {
            Re = std::stod(argv[++i]);
        } else if (arg == "--trt") {
            use_trt = true;
        } else if (arg == "--magic" && i + 1 < argc) {
            magic_param = std::stod(argv[++i]);
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_usage();
            return 1;
        }
    }

    // --- Create and Run Test Case ---
    try {
        // Create test case
        std::unique_ptr<TestCase> test_case;
        try {
            test_case = createTestCase(test_case_name);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            print_usage();
            return 1;
        }

        std::cout << "Running test case: " << test_case->getName()
                  << " - " << test_case->getDescription() << std::endl;

        // Set up simulator with test case parameters
        std::unique_ptr<Simulator> sim = test_case->setup(nx, ny, Re, use_trt, magic_param);

        // Run simulation
        sim->run();

    } catch (const std::exception& e) {
        std::cerr << "Error: Simulation failed with exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Error: Simulation failed with unknown exception." << std::endl;
        return 1;
    }

    std::cout << "Simulation completed successfully." << std::endl;
    return 0;
}