# Simulate zero coherence trials for selected sigmas in order to illustrate
# the different regimes (5000 trials for each sigma).
# These simulations are used for computing the PKs.

COH=0
SIGMA=40
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 1000 &
./start.sh 1 $COH $SIGMA 1001 2000 &
./start.sh 1 $COH $SIGMA 2001 3000 &
./start.sh 1 $COH $SIGMA 3001 4000 &
./start.sh 1 $COH $SIGMA 4001 5000 &
cd ..

COH=0
SIGMA=80
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 1000 &
./start.sh 1 $COH $SIGMA 1001 2000 &
./start.sh 1 $COH $SIGMA 2001 3000 &
./start.sh 1 $COH $SIGMA 3001 4000 &
./start.sh 1 $COH $SIGMA 4001 5000 &
cd ..

COH=0
SIGMA=140
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_all_to_all.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 1000 &
./start.sh 1 $COH $SIGMA 1001 2000 &
./start.sh 1 $COH $SIGMA 2001 3000 &
./start.sh 1 $COH $SIGMA 3001 4000 &
./start.sh 1 $COH $SIGMA 4001 5000 &
cd ..

