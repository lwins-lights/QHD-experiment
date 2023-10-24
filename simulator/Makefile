all: pseudospec nagd sgd

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

sgd.o: sgd.cpp config.hpp
	$(CXX) $(CXXFLAGS) -c $^ -fopenmp

sgd: sgd.o potential.o
	$(CXX) $(CXXFLAGS) -o sgd $^ -fopenmp -lm -lcnpy -lz -ltbb

qipopt_F.o: qipopt_F.hpp qipopt_F.cpp
	$(CXX) $(CXXFLAGS) -c $^

qipopt.o: qipopt.cpp config.hpp
	$(CXX) $(CXXFLAGS) -c $^ -fopenmp

qipopt: qipopt.o potential.o
	$(CXX) $(CXXFLAGS) -o qipopt $^ -fopenmp -lfftw3_omp -lfftw3 -lm -lcnpy -lz -ltbb

clean:
	rm -f pseudospec pseudospec.o nagd nagd.o sgd sgd.o potential.hpp.gch potential.o config.hpp.gch qipopt qipopt.o