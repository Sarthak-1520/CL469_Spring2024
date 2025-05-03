# Compiler
CXX = g++
# Compiler flags: -O3 for optimization, -Wall for warnings
CXXFLAGS = -std=c++17 -O3 -Wall -Iinclude
# Linker flags
LDFLAGS = -lm # Link math library

# Source files directory
SRCDIR = src
# Include files directory
INCDIR = include
# Object files directory
OBJDIR = obj
# Executable name
EXECUTABLE = entropic_lbm
# Results directory
RESDIR = results
# Figures directory
FIGDIR = figures

# Find all .cpp files in src directory
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
# Create object file names by replacing src/%.cpp with obj/%.o
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Default target
all: $(EXECUTABLE)

# Link the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
	@echo "Executable '$@' created."

# Compile source files into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create object directory if it doesn't exist
$(OBJDIR):
	@mkdir -p $(OBJDIR)

# Create results directory if it doesn't exist
$(RESDIR):
	@mkdir -p $(RESDIR)

# Create figures directory if it doesn't exist
$(FIGDIR):
	@mkdir -p $(FIGDIR)

# Target to run the simulation and plot convergence
run: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running simulation..."
	./$(EXECUTABLE)
	@echo "Simulation finished. Visualizing results with Python..."
	@python visualize_simple.py --test ellipse_flow --model entropic
	@echo "Python visualization script finished."

# Target to run the simulation and generate data files only (no visualization)
run-data-only: $(EXECUTABLE) $(RESDIR)
	@echo "Running simulation (data only)..."
	./$(EXECUTABLE)
	@echo "Simulation finished. Data files generated in $(RESDIR)/"

# Test targets for different test cases
test-ellipse: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running elliptical flow test case with entropic collision..."
	./$(EXECUTABLE) --test ellipse_flow
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test ellipse_flow --model entropic
	@echo "Python visualization script finished."
	@echo "Running elliptical flow test case with TRT collision..."
	./$(EXECUTABLE) --test ellipse_flow --trt
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test ellipse_flow --model trt
	@echo "Python visualization script finished."

test-cavity: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running lid-driven cavity test case with entropic collision..."
	./$(EXECUTABLE) --test lid_driven_cavity
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test lid_driven_cavity --model entropic
	@echo "Python visualization script finished."
	@echo "Running lid-driven cavity test case with TRT collision..."
	./$(EXECUTABLE) --test lid_driven_cavity --trt
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test lid_driven_cavity --model trt
	@echo "Python visualization script finished."

test-channel: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running channel flow with obstacle test case with entropic collision..."
	./$(EXECUTABLE) --test channel_obstacle
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test channel_obstacle --model entropic
	@echo "Python visualization script finished."
	@echo "Running channel flow with obstacle test case with TRT collision..."
	./$(EXECUTABLE) --test channel_obstacle --trt
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test channel_obstacle --model trt
	@echo "Python visualization script finished."

test-step: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running backward-facing step test case with entropic collision..."
	./$(EXECUTABLE) --test backward_facing_step
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test backward_facing_step --model entropic
	@echo "Python visualization script finished."
	@echo "Running backward-facing step test case with TRT collision..."
	./$(EXECUTABLE) --test backward_facing_step --trt
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test backward_facing_step --model trt
	@echo "Python visualization script finished."

test-vortex: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running Taylor-Green vortex test case with entropic collision..."
	./$(EXECUTABLE) --test taylor_green_vortex
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test taylor_green_vortex --model entropic
	@echo "Python visualization script finished."
	@echo "Running Taylor-Green vortex test case with TRT collision..."
	./$(EXECUTABLE) --test taylor_green_vortex --trt
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test taylor_green_vortex --model trt
	@echo "Python visualization script finished."

test-poiseuille: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running Poiseuille flow test case with entropic collision..."
	./$(EXECUTABLE) --test poiseuille_flow
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test poiseuille_flow --model entropic
	@echo "Python visualization script finished."
	@echo "Running Poiseuille flow test case with TRT collision..."
	./$(EXECUTABLE) --test poiseuille_flow --trt
	@echo "Test finished. Visualizing results with Python..."
	@python visualize_simple.py --type all --test poiseuille_flow --model trt
	@echo "Python visualization script finished."

# Compare entropic and TRT collision models
compare-models: $(EXECUTABLE) $(RESDIR) $(FIGDIR)
	@echo "Running elliptical flow with entropic collision..."
	./$(EXECUTABLE) --test ellipse_flow
	@echo "Renaming results directory for entropic model..."
	@mv $(RESDIR) $(RESDIR)_entropic
	@mkdir -p $(RESDIR)
	@echo "Running elliptical flow with TRT collision..."
	./$(EXECUTABLE) --test ellipse_flow --trt
	@echo "Renaming results directory for TRT model..."
	@mv $(RESDIR) $(RESDIR)_trt
	@echo "Comparison complete. Results in $(RESDIR)_entropic and $(RESDIR)_trt."

# Clean up object files and the executable
clean:
	@echo "Cleaning up..."
	@rm -rf $(OBJDIR) $(EXECUTABLE) $(RESDIR)/* $(FIGDIR)/* core.* *.dSYM

# Clean up results only
clean-results:
	@echo "Cleaning up results..."
	@rm -rf $(RESDIR)/* $(FIGDIR)/*

# Generate .dat files only (no visualization)
data-ellipse: $(EXECUTABLE) $(RESDIR)
	@echo "Running elliptical flow test case with entropic collision (data only)..."
	./$(EXECUTABLE) --test ellipse_flow
	@echo "Data files generated in $(RESDIR)/ellipse_flow_entropic/"
	@echo "Running elliptical flow test case with TRT collision (data only)..."
	./$(EXECUTABLE) --test ellipse_flow --trt
	@echo "Data files generated in $(RESDIR)/ellipse_flow_trt/"

data-cavity: $(EXECUTABLE) $(RESDIR)
	@echo "Running lid-driven cavity test case with entropic collision (data only)..."
	./$(EXECUTABLE) --test lid_driven_cavity
	@echo "Data files generated in $(RESDIR)/lid_driven_cavity_entropic/"
	@echo "Running lid-driven cavity test case with TRT collision (data only)..."
	./$(EXECUTABLE) --test lid_driven_cavity --trt
	@echo "Data files generated in $(RESDIR)/lid_driven_cavity_trt/"

data-channel: $(EXECUTABLE) $(RESDIR)
	@echo "Running channel flow with obstacle test case with entropic collision (data only)..."
	./$(EXECUTABLE) --test channel_obstacle
	@echo "Data files generated in $(RESDIR)/channel_obstacle_entropic/"
	@echo "Running channel flow with obstacle test case with TRT collision (data only)..."
	./$(EXECUTABLE) --test channel_obstacle --trt
	@echo "Data files generated in $(RESDIR)/channel_obstacle_trt/"

data-step: $(EXECUTABLE) $(RESDIR)
	@echo "Running backward-facing step test case with entropic collision (data only)..."
	./$(EXECUTABLE) --test backward_facing_step
	@echo "Data files generated in $(RESDIR)/backward_facing_step_entropic/"
	@echo "Running backward-facing step test case with TRT collision (data only)..."
	./$(EXECUTABLE) --test backward_facing_step --trt
	@echo "Data files generated in $(RESDIR)/backward_facing_step_trt/"

data-vortex: $(EXECUTABLE) $(RESDIR)
	@echo "Running Taylor-Green vortex test case with entropic collision (data only)..."
	./$(EXECUTABLE) --test taylor_green_vortex
	@echo "Data files generated in $(RESDIR)/taylor_green_vortex_entropic/"
	@echo "Running Taylor-Green vortex test case with TRT collision (data only)..."
	./$(EXECUTABLE) --test taylor_green_vortex --trt
	@echo "Data files generated in $(RESDIR)/taylor_green_vortex_trt/"

data-poiseuille: $(EXECUTABLE) $(RESDIR)
	@echo "Running Poiseuille flow test case with entropic collision (data only)..."
	./$(EXECUTABLE) --test poiseuille_flow
	@echo "Data files generated in $(RESDIR)/poiseuille_flow_entropic/"
	@echo "Running Poiseuille flow test case with TRT collision (data only)..."
	./$(EXECUTABLE) --test poiseuille_flow --trt
	@echo "Data files generated in $(RESDIR)/poiseuille_flow_trt/"

# Generate all data files
data-all: data-ellipse data-cavity data-channel data-step data-vortex data-poiseuille
	@echo "All data files generated in $(RESDIR)/"

# Phony targets
.PHONY: all clean clean-results run run-data-only test-ellipse test-cavity test-channel test-step test-vortex test-poiseuille compare-models data-ellipse data-cavity data-channel data-step data-vortex data-poiseuille data-all