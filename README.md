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
First, define the potential function in a `.cpp` file (refer to `func/l2.cpp` for example). The following shell script will run QHD and NAGD to find global minima of the potential function:
`bash ./run_all.sh <num_cells> <T> <dt> <func_cpp_path>`
-   `num_cells` dictates the number of cells in each dimension due to spatial discretization in QHD. 
-   `T` is the evolution time.
-   `dt` is the time step for each iteration due to time discretization. 
- `func_cpp_path` gives the path of the `.cpp` file for the potential function.
- 
Example:
`bash ./run_all.sh 256 10 0.001 func/ackley2.cpp`

To visualize the result, simply use the Python script:
`python3 plot_result.py`

## Credit 
Some code is adapted form https://github.com/jiaqileng/quantum-hamiltonian-descent/.