# bash run_all.sh <fcp> <len> <tot_time> <dt> <par> <tot_steps> <learning_rate> <sample_number> <noise_level> <L>
fcp=$1;
len=$2;
tot=$3;
dt=$4;
par=$5;
ts=$6;
lr=$7;
sn=$8;
nl=$9;
L=${10};

mkdir -p result
cp $fcp simulator/potential.cpp
cd simulator
make pseudospec subgrad lfmsgd
./pseudospec $len $tot $dt $par $L
./subgrad $ts $lr $par $sn $L
./lfmsgd $ts $nl $par $sn $L
make clean
cd ..

