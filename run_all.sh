# bash run_all.sh <len> <tot_time> <dt> <par> <noise_level> <func_cpp_path>
len=$1;
tot=$2;
dt=$3;
par=$4;
nl=$5;
fcp=$6;

mkdir -p result
cp $fcp simulator/potential.cpp
cd simulator
make
./pseudospec $len $tot $dt $par
./nagd $len $tot $dt $par
./sgd $len $tot $dt $par $nl
make clean
cd ..

