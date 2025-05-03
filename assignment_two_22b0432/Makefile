# Makefile for CL469 Assignment Two

# Compiler
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Iinclude

# Directories
SRCDIR = src
INCLUDEDIR = include
BUILDDIR = build
DATADIR = data
VISDIR = visualizations
FIGDIR = figures
REPORT = report

# Source files and object files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))

# Executable name
EXECUTABLE = $(BUILDDIR)/lbm_poiseuille

# Default target
all: $(EXECUTABLE)

# Link object files to create executable
$(EXECUTABLE): $(OBJECTS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compile source files into object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build and data directories if they don't exist
$(BUILDDIR) $(DATADIR):
	mkdir -p $@

# Target to run the simulation
run: $(EXECUTABLE) $(DATADIR)
	./$(EXECUTABLE)

# Target to generate visualizations (Placeholder)
visualize:
	# Add commands to run python scripts in $(VISDIR)
	@echo "Generating visualizations..."
	python $(VISDIR)/plot_velocity.py
	python $(VISDIR)/plot_forces.py
	python $(VISDIR)/plot_convergence.py
	# Run the report figure generation script - requires data from multiple runs!
	# Comment out if data subdirectories are not set up
	python $(VISDIR)/plot_report_figures.py
	# ... add other visualization scripts

# Target to compile the LaTeX report
report: $(REPORT).tex $(FIGDIR)/*
	pdflatex $(REPORT).tex
	pdflatex $(REPORT).tex # Run twice for references/toc
	@echo "Report generated: $(REPORT).pdf"

# Target to clean build files, data, figures, and executable
clean:
	rm -rf $(BUILDDIR) $(EXECUTABLE) $(DATADIR)/* $(FIGDIR)/* $(REPORT).aux $(REPORT).log $(REPORT).pdf $(REPORT).toc

# Phony targets
.PHONY: all run visualize report clean 