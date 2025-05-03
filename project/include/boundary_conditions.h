#pragma once
#include "node_data.h"
#include "lattice.h"
#include <vector>

// Forward declaration
class Simulator;
enum class GeometryType {
    ELLIPSE,               // Flow around an elliptical obstacle
    LID_DRIVEN_CAVITY,     // Lid-driven cavity flow
    CHANNEL_OBSTACLE,      // Channel flow with a circular obstacle
    BACKWARD_FACING_STEP,  // Backward-facing step flow
    TAYLOR_GREEN_VORTEX,   // Taylor-Green vortex decay
    POISEUILLE_FLOW        // Poiseuille flow in a channel
};

class BoundaryConditions {
public:
    /**
     * @brief Apply all boundary conditions to the grid
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     * @param characteristic_velocity The characteristic velocity (inlet or lid)
     * @param geometry_type The type of geometry being simulated
     */
    static void apply_all(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice,
        Real characteristic_velocity,
        GeometryType geometry_type = GeometryType::ELLIPSE
    );

    /**
     * @brief Apply bounce-back on all wall nodes
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     */
    static void apply_bounce_back_walls(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice
    );

    /**
     * @brief Apply inlet boundary condition
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     * @param inlet_velocity The inlet velocity
     */
    static void apply_inlet_boundary(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice,
        Real inlet_velocity
    );

    /**
     * @brief Apply outlet boundary condition
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     */
    static void apply_outlet_boundary(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice
    );

    /**
     * @brief Apply lid boundary condition
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     * @param lid_velocity The lid velocity
     * @param nx Grid size in x direction
     * @param ny Grid size in y direction
     */
    static void apply_lid_boundary(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice,
        Real lid_velocity,
        int nx, int ny
    );

    /**
     * @brief Apply partial inlet boundary condition
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     * @param inlet_velocity The inlet velocity
     * @param nx Grid size in x direction
     * @param ny Grid size in y direction
     */
    static void apply_partial_inlet_boundary(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice,
        Real inlet_velocity,
        int nx, int ny
    );

    /**
     * @brief Apply pressure boundary conditions
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     * @param nx Grid size in x direction
     * @param ny Grid size in y direction
     */
    static void apply_pressure_boundaries(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice,
        int nx, int ny
    );

    /**
     * @brief Simple bounce-back for solid walls
     *
     * @param wall_node The wall node (assumes node.is_fluid == false)
     * @param lattice The lattice definition
     * @param fluid_neighbors The neighboring fluid nodes
     */
    static void bounce_back(
        NodeData& wall_node,
        const Lattice& lattice,
        const std::vector<NodeData*>& fluid_neighbors
    );

    /**
     * @brief Fixed velocity boundary (Zou/He style)
     *
     * @param boundary_node The boundary node
     * @param lattice The lattice definition
     * @param target_u The target velocity
     * @param target_rho The target density (negative means keep node density)
     */
    static void fixed_velocity(
        NodeData& boundary_node,
        const Lattice& lattice,
        const Vector2D& target_u,
        Real target_rho = -1.0
    );

    /**
     * @brief Periodic boundary condition
     *
     * @param grid The simulation grid
     * @param lattice The lattice definition
     * @param nx Grid size in x direction
     * @param ny Grid size in y direction
     */
    static void apply_periodic_boundaries(
        std::vector<std::vector<NodeData>>& grid,
        const Lattice& lattice,
        int nx, int ny
    );

    /**
     * @brief Pressure boundary condition (density-based)
     *
     * @param boundary_node The boundary node
     * @param lattice The lattice definition
     * @param target_rho The target density
     * @param is_inlet Whether this is an inlet (true) or outlet (false)
     */
    static void pressure_boundary(
        NodeData& boundary_node,
        const Lattice& lattice,
        Real target_rho,
        bool is_inlet
    );
};