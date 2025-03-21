# QHD-experiment

This is the repository for [Quantum Hamiltonian Descent for Non-smooth Optimization (arXiv:2503.15878)](https://arxiv.org/abs/2503.15878).

## Dependencies
```
sudo apt install make g++ git libtbb-dev python3-pip
pip install matplotlib scikit-optimize tqdm
```

### C++ libraries
It is required to install [fftw](https://www.fftw.org) (with OpenMP enabled) and [cnpy](https://github.com/rogersce/cnpy).
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

Make sure you have a enough stack size by setting `ulimit -s unlimited`.

### Benchmarking QHD against Subgrad and LFMSGD

Simply run `python3 script/new_skopt_find_best_L.py` to obtain the results that will be aggregated in `result/data.csv`. Note that unlike in the paper, all three algorithms QHD, Subgrad and LFMSGD are parameterized by $L$ in the code: it is easy to see that optimizing $L$ is equivalent to optimizing $\eta$ and $\sigma$ in the paper.

### Visualizing the landscape of a 2D test function

This can be done by running `python3 script/show_discretized_potential_2d.py`. The resulting figure is saved as `result/landscape.png`.

### Visualizing distribution at intermediate iterations of QHD, Subgrad and LFMSGD

This is handled by `python3 script/visualize_dist_2d.py` and the output images are saved in `result/visualize_dist_2d/`.

## Credit 
Some code is adapted from https://github.com/jiaqileng/quantum-hamiltonian-descent/.
