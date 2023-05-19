all: pseudospec nagd

CXX := g++
CXXFLAGS := -O3 --std=c++17

potential.o: potential.hpp potential.cpp
	$(CXX) $(CXXFLAGS) -c $^

pseudospec.o: pseudospec.cpp config.hpp
	$(CXX) $(CXXFLAGS) -c $^ -fopenmp

pseudospec: pseudospec.o potential.o
	$(CXX) $(CXXFLAGS) -o pseudospec $^ -fopenmp -lfftw3_omp -lfftw3 -lm -lcnpy -lz -ltbb

nagd.o: nagd.cpp config.hpp
	$(CXX) $(CXXFLAGS) -c $^ -fopenmp

nagd: nagd.o potential.o
	$(CXX) $(CXXFLAGS) -o nagd $^ -fopenmp -lm -lcnpy -lz -ltbb

clean:
	rm -f pseudospec pseudospec.o nagd nagd.o potential.hpp.gch potential.o config.hpp.gch