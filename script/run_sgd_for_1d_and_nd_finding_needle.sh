#!/bin/bash

BLUE='\033[0;34m'
NC='\033[0m' # No Color

# run_all.sh <len> <tot_time> <dt> <par> <noise_level> <func_cpp_path>

# run <func_name> <len> <tot_time> <dt> <par> <noise_level> <max_sample_number>
run() {
    echo -e "${BLUE}Problem $1 (noise_level=$6)${NC}"
    mkdir -p result
    cp "func/$1.cpp" simulator/potential.cpp
    cd simulator
    make
    ./sgd $2 $3 $4 $5 $6 $7
    cd .. 
    python3 plot_result.py -s -S
    cp result/result.png "result/1d_and_nd_find_needle/${1}(${6}).png"
}

mkdir -p result/1d_and_nd_find_needle

run nd_find_needle 4096 50 0.001 1 0.96 10000
run nd_find_needle 4096 50 0.001 1 1.12 10000
run nd_find_needle 4096 50 0.001 1 1.28 10000
run nd_find_needle 4096 50 0.001 1 1.44 10000
run 1d_find_needle 4096 50 0.001 1 0.48 10000
run 1d_find_needle 4096 50 0.001 1 0.56 10000
run 1d_find_needle 4096 50 0.001 1 0.64 10000
run 1d_find_needle 4096 50 0.001 1 0.72 10000
