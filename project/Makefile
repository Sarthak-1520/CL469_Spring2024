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

# Target to run the simulation and plot convergence
run: $(EXECUTABLE) $(RESDIR)
	@echo "Running simulation..."
	./$(EXECUTABLE)
	@echo "Simulation finished. Visualizing results with Python..."
	@python visualize_results.py # Changed from gnuplot
	@echo "Python visualization script finished."

# Clean up object files and the executable
clean:
	@echo "Cleaning up..."
	@rm -rf $(OBJDIR) $(EXECUTABLE) $(RESDIR)/* core.* *.dSYM

# Phony targets
.PHONY: all clean run