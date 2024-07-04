# bash run_all.sh <fcp> <len> <tot_time> <dt> <par> <tot_steps> <learning_rate> <sample_number>
fcp=$1;
len=$2;
tot=$3;
dt=$4;
par=$5;
ts=$6;
lr=$7;
sn=$8;

mkdir -p result
cp $fcp simulator/potential.cpp
cd simulator
make pseudospec subgrad
./pseudospec $len $tot $dt $par
./subgrad $ts $lr $par $sn 
make clean
cd ..

