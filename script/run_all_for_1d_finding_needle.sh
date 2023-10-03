#!/bin/bash

BLUE='\033[0;34m'
NC='\033[0m' # No Color

# run_all.sh <len> <tot_time> <dt> <par> <noise_level> <func_cpp_path>

# run <func_name> <len> <tot_time> <dt> <par> <noise_level>
run() {
    echo -e "${BLUE}Problem $1 (noise_level=$6)${NC}"
    bash run_all.sh $2 $3 $4 $5 $6 "func/$1.cpp"
    python3 plot_result.py -s
    cp result/result.png "result/1d_find_needle/${1}(${6}).png"
}

mkdir -p result/1d_find_needle

#run 1d_find_needle 4096 200 0.001 1 0.01
#run 1d_find_needle 4096 200 0.001 1 0.02
#run 1d_find_needle 4096 200 0.001 1 0.04
#run 1d_find_needle 4096 200 0.001 1 0.08
#run 1d_find_needle 4096 200 0.001 1 0.16
#run 1d_find_needle 4096 200 0.001 1 0.32
#run 1d_find_needle 4096 200 0.001 1 0.64
#run 1d_find_needle 4096 200 0.001 1 1.28

#run 1d_find_needle 4096 200 0.001 1 0.40
#run 1d_find_needle 4096 200 0.001 1 0.48
#run 1d_find_needle 4096 200 0.001 1 0.56
#run 1d_find_needle 4096 200 0.001 1 0.72
#run 1d_find_needle 4096 200 0.001 1 0.80
#run 1d_find_needle 4096 200 0.001 1 0.88
#run 1d_find_needle 4096 200 0.001 1 0.96
#run 1d_find_needle 4096 200 0.001 1 1.04
#run 1d_find_needle 4096 200 0.001 1 1.12
#run 1d_find_needle 4096 200 0.001 1 1.20

run 1d_find_needle 4096 1000 0.001 1 0.46
run 1d_find_needle 4096 1000 0.001 1 0.47
run 1d_find_needle 4096 1000 0.001 1 0.48
run 1d_find_needle 4096 1000 0.001 1 0.49
run 1d_find_needle 4096 1000 0.001 1 0.50