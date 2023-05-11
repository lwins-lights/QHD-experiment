# QHD-experiment

## Dependencies
```
sudo apt install make g++ python3-pip
pip install matplotlib
```

### C++ libraries
It is required to install [`fftw`](https://www.fftw.org) (with OpenMP enabled) and [`cnpy`](https://github.com/rogersce/cnpy).
For `fftw`:

```
# download fftw
wget https://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvf fftw-3.3.10.tar.gz
cd fftw-3.3.10

# install fftw with OpenMP enabled
./configure --enable-openmp
make
sudo make install
```
For `cnpy`:
```
# install dependencies
sudo apt install cmake zlib1g-dev

# download cnpy
git clone https://github.com/rogersce/cnpy.git
cd cnpy

# install cnpy
mkdir build
cd build
cmake ..
make
sudo make install

# configure for the dymanic library just installed
sudo ldconfig
```

## Usage
First, define the potential function in a `.cpp` file (refer to `func/l2.cpp` for example). The potential will be defined on [`-L`,`L`)^`dim` where `L` and `dim` specified in the corresponding `.cpp` file.
The following shell script will run solvers (QHD and NAGD) to find the global minimum of the potential function:
`bash ./run_all.sh <num_cells> <T> <dt> <par> <func_cpp_path>`
-   `num_cells` dictates the number of cells in each dimension due to spatial discretization in QHD. 
-   `T` is the evolution time.
-   `dt` is the time step for each iteration due to time discretization. 
-   `func_cpp_path` gives the path of the `.cpp` file for the potential function.
-   `par` tells the parallelism number. For example, `par = 2` will make all solvers run twice, and the final result for each solver will be the smaller one of those of two runs.

Example:
`bash ./run_all.sh 256 10 0.001 1 func/ackley2.cpp`

To visualize the result, simply use the Python script:
`python3 plot_result.py`

## Possible Issues
### Segmentation fault
This is likely due to an insufficient stack size which is set to 8MB by default. Use, say, `ulimit -s 262144` to increase it to 256MB.

## Credit 
Some code is adapted from https://github.com/jiaqileng/quantum-hamiltonian-descent/.