all: pseudospec subgrad lfmsgd

CXX := g++
CXXFLAGS := -O3 --std=c++17

#potential.o: potential.cpp potential_util.cpp potential.hpp
#	$(CXX) $(CXXFLAGS) -c $^

#pseudospec.o: pseudospec.cpp config.hpp
#	$(CXX) $(CXXFLAGS) -c $^ -fopenmp

#pseudospec: pseudospec.o potential.o
#	$(CXX) $(CXXFLAGS) -o pseudospec $^ -fopenmp -lfftw3_omp -lfftw3 -lm -lcnpy -lz -ltbb

pseudospec: pseudospec.cpp config.hpp potential.cpp potential_util.cpp potential.hpp
	$(CXX) $(CXXFLAGS) -o pseudospec $^ -fopenmp -lfftw3_omp -lfftw3 -lm -lcnpy -lz -ltbb

subgrad: subgrad.cpp config.hpp potential.cpp potential_util.cpp potential.hpp
	$(CXX) $(CXXFLAGS) -o subgrad $^ -fopenmp -lm -lcnpy -lz -ltbb

lfmsgd: lfmsgd.cpp config.hpp potential.cpp potential_util.cpp potential.hpp
	$(CXX) $(CXXFLAGS) -o lfmsgd $^ -fopenmp -lm -lcnpy -lz -ltbb

clean:
	rm -f pseudospec subgrad lfmsgd