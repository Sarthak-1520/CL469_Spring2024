CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11
LATEX = pdflatex
PYTHON = python3

all: simulation visualize report

simulation: poiseuille_flow.cpp
	mkdir -p output
	$(CXX) $(CXXFLAGS) -o poiseuille_flow poiseuille_flow.cpp
	./poiseuille_flow

visualize: simulation plot_results.py
	$(PYTHON) plot_results.py

report: report.tex
	$(LATEX) report.tex
	$(LATEX) report.tex
	
animate: visualize
	$(PYTHON) animate_profiles.py

clean:
	rm -f poiseuille_flow *.o *.aux *.log *.out
	rm -rf output/
	rm -rf figures/
	rm -f animate_profiles.py
	rm -f velocity_evolution.mp4

.PHONY: all clean simulation report visualize animate 