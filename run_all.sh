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

# example command:
# python script/find_best_L.py --fpath func/nonsmooth/.cpp 
# python script/plot_all.py --fpath func/nonsmooth/.cpp 
# python script/run_qhd.py --fpath func/nonsmooth/.cpp --Llist 0.1,0.2,0.4,1,2,4,10,20,50,100
# python script/run_qhd.py --fpath func/nonsmooth/.cpp --Llist 1,2,4,10,25,50,100,250,500,1000,1500,2000
# python script/run_gradtest.py --fpath func/nonsmooth/.cpp