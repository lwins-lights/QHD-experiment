#!/bin/bash

# run_all.sh <len> <tot_time> <dt> <par> <func_cpp_path>

# run <func_name> <len> <tot_time> <dt> <par>
run() {
    bash run_all.sh $2 $3 $4 $5 "func/floudas2013handbook/$1.cpp"
    python3 plot_result.py -s
    cp result/result.png result/floudas2013handbook/$1.png
}

mkdir -p result/floudas2013handbook

run p4_5 256 10 0.0001 1
run p4_6 65536 10 0.0001 1

