# bash run_all.sh <len> <tot_time> <dt> <func_cpp_path>
len=$1;
tot=$2;
dt=$3;
fcp=$4;

mkdir -p result
cp $fcp simulator/potential.cpp
cd simulator
make
./pseudospec $len $tot $dt
./nagd $len $tot $dt
make clean
cd ..

