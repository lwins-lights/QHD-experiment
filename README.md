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

### Optimizing for one function
First, define the potential function in a `.cpp` file (refer to `func/nonsmooth/l2.cpp` for example). The potential will be defined on [`-L`,`L`)^`dim` where `L` and `dim` specified in the corresponding `.cpp` file.
The following shell script will run solvers (QHD and NAGD) to find the global minimum of the potential function:
`bash ./run_all.sh <func_cpp_path> <num_cells> <T> <dt> <par> <tot_steps> <learning_rate> <sample_number>`
-   `func_cpp_path` gives the path of the `.cpp` file for the potential function.
-   `num_cells` dictates the number of cells in each dimension due to spatial discretization in QHD. 
-   `T` is the evolution time.
-   `dt` is the time step for each iteration due to time discretization. 
-   `par` tells the parallelism number. For example, `par = 2` will make all solvers run twice, and the final result for each solver will be the smaller one of those of two runs.
-   `tot_steps` is the total number of steps for the subgradient method solver.
-   `learning_rate` controls the learning rate of the subgradient method. *The recommended value is $L/G$ where $L$ is `L` and $G$ is the Lipschitz constant of the potential function.*
-   `sample_number` is the number of samples used in subgradient method solver for deriving probability distributions.

Example:
`./run_all.sh func/nonsmooth/l2.cpp 256 10 0.001 1 10000 1 10000`

To visualize the result, simply use the Python script:
`python3 plot_result.py`

### Preset script
Use `bash script/run_all_in_floudas2013handbook.sh`, which in turn calls `run_all.sh`, to run the preset experiment for selected functions from *Handbook of Test Problems in Local and Global Optimization*. The result will be saved in `result/floudas2013handbook/`.

## Possible Issues
### Segmentation fault
This is likely due to an insufficient stack size which is set to 8MB by default. Use, say, `ulimit -s 262144` to increase it to 256MB.

## Credit 
Some code is adapted from https://github.com/jiaqileng/quantum-hamiltonian-descent/.
