#!/bin/bash

BLUE='\033[0;34m'
NC='\033[0m' # No Color

# run_all.sh <len> <tot_time> <dt> <par> <func_cpp_path>

# run <func_name> <len> <tot_time> <dt> <par>
run() {
    echo -e "${BLUE}Problem $1${NC}"
    bash run_all.sh $2 $3 $4 $5 "func/floudas2013handbook/$1.cpp"
    python3 plot_result.py -s
    cp result/result.png result/floudas2013handbook/$1.png
}

mkdir -p result/floudas2013handbook

# 1D funcs 
#run p4_3 65536 10 0.0001 1
#run p4_4 65536 10 0.0001 1
#run p4_6 65536 10 0.0001 1
#run p4_7 65536 10 0.0001 1

# 2D funcs
#run p4_5 256 10 0.0001 1
#run p4_9 256 10 0.0001 1

# 3D funcs
#run p7_6 32 10 0.0001 1

# 5D funcs
run p2_1 8 10 0.0001 1