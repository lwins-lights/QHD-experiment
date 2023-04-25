all: pseudospec

CXX := g++
CXXFLAGS := -O3 --std=c++17

potential.o: potential.hpp potential.cpp
	$(CXX) $(CXXFLAGS) -c $^

pseudospec.o: pseudospec.cpp
	$(CXX) $(CXXFLAGS) -c $^ -fopenmp

pseudospec: pseudospec.o potential.o
	$(CXX) $(CXXFLAGS) -o pseudospec $^ -fopenmp -lfftw3_omp -lfftw3 -lm

clean:
	rm -f pseudospec pseudospec.o potential.hpp.gch potential.o 