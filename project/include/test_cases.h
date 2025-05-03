#pragma once
#include "common_types.h"
#include "simulator.h"
#include <string>
#include <functional>
#include <memory>

/**
 * @brief Abstract base class for test cases
 * 
 * This class defines the interface for all test cases.
 * Each test case must implement the setup method to configure
 * the simulation parameters and boundary conditions.
 */
class TestCase {
public:
    virtual ~TestCase() = default;
    
    /**
     * @brief Set up the test case
     * 
     * @param nx Grid size in x direction
     * @param ny Grid size in y direction
     * @param Re Reynolds number
     * @param use_trt Whether to use TRT collision model
     * @param magic_param Magic parameter for TRT model
     * @return std::unique_ptr<Simulator> Configured simulator
     */
    virtual std::unique_ptr<Simulator> setup(
        int nx, int ny, Real Re, bool use_trt = false, Real magic_param = 0.25) = 0;
    
    /**
     * @brief Get the name of the test case
     * 
     * @return std::string Name of the test case
     */
    virtual std::string getName() const = 0;
    
    /**
     * @brief Get the description of the test case
     * 
     * @return std::string Description of the test case
     */
    virtual std::string getDescription() const = 0;
};

/**
 * @brief Flow around an elliptical obstacle
 */
class EllipseFlowTestCase : public TestCase {
public:
    std::unique_ptr<Simulator> setup(
        int nx, int ny, Real Re, bool use_trt = false, Real magic_param = 0.25) override;
    
    std::string getName() const override {
        return "ellipse_flow";
    }
    
    std::string getDescription() const override {
        return "Flow around an elliptical obstacle";
    }
};

/**
 * @brief Lid-driven cavity flow
 */
class LidDrivenCavityTestCase : public TestCase {
public:
    std::unique_ptr<Simulator> setup(
        int nx, int ny, Real Re, bool use_trt = false, Real magic_param = 0.25) override;
    
    std::string getName() const override {
        return "lid_driven_cavity";
    }
    
    std::string getDescription() const override {
        return "Lid-driven cavity flow";
    }
};

/**
 * @brief Channel flow with a circular obstacle
 */
class ChannelObstacleTestCase : public TestCase {
public:
    std::unique_ptr<Simulator> setup(
        int nx, int ny, Real Re, bool use_trt = false, Real magic_param = 0.25) override;
    
    std::string getName() const override {
        return "channel_obstacle";
    }
    
    std::string getDescription() const override {
        return "Channel flow with a circular obstacle";
    }
};

/**
 * @brief Backward-facing step flow
 */
class BackwardFacingStepTestCase : public TestCase {
public:
    std::unique_ptr<Simulator> setup(
        int nx, int ny, Real Re, bool use_trt = false, Real magic_param = 0.25) override;
    
    std::string getName() const override {
        return "backward_facing_step";
    }
    
    std::string getDescription() const override {
        return "Backward-facing step flow";
    }
};

/**
 * @brief Taylor-Green vortex decay
 */
class TaylorGreenVortexTestCase : public TestCase {
public:
    std::unique_ptr<Simulator> setup(
        int nx, int ny, Real Re, bool use_trt = false, Real magic_param = 0.25) override;
    
    std::string getName() const override {
        return "taylor_green_vortex";
    }
    
    std::string getDescription() const override {
        return "Taylor-Green vortex decay";
    }
};

/**
 * @brief Poiseuille flow in a channel
 */
class PoiseuilleFlowTestCase : public TestCase {
public:
    std::unique_ptr<Simulator> setup(
        int nx, int ny, Real Re, bool use_trt = false, Real magic_param = 0.25) override;
    
    std::string getName() const override {
        return "poiseuille_flow";
    }
    
    std::string getDescription() const override {
        return "Poiseuille flow in a channel";
    }
};

/**
 * @brief Factory function to create a test case by name
 * 
 * @param name Name of the test case
 * @return std::unique_ptr<TestCase> Test case instance
 */
std::unique_ptr<TestCase> createTestCase(const std::string& name);

/**
 * @brief Get a list of all available test cases
 * 
 * @return std::vector<std::string> List of test case names
 */
std::vector<std::string> getAvailableTestCases();
