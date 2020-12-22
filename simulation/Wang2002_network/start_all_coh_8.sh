# Simulate 5000 trials for a fixed coherence and varying magnitude of stimulus fluctuations
# sigma in order to measure the accuracy as a function of sigma.

COH=8
SIGMA=0
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=8
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=16
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=24
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=32
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=40
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=48
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=56
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=64
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=72
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=80
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=88
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=96
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=104
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..


COH=8
SIGMA=112
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=120
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=128
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=136
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=144
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=152
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=8
SIGMA=160
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..
