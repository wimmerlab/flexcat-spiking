# Simulate zero coherence trials for selected sigmas in order to illustrate
# the different regimes (5000 trials for each sigma).
# These simulations are used for computing the PKs.

COH=0
SIGMA=8
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=0
SIGMA=20
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=0
SIGMA=36
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

