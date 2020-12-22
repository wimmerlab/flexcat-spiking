# Simulate 5000 trials for a fixed coherence and varying magnitude of stimulus fluctuations
# sigma in order to measure the accuracy as a function of sigma.

COH=1
SIGMA=0
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=4
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=8
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=12
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=16
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=20
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=24
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=28
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=32
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=36
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=40
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=44
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=48
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

COH=1
SIGMA=52
MYDIR=coh_${COH}_sigma_${SIGMA}
mkdir $MYDIR
cp decision_network_sparse.py $MYDIR/worker1.py
cp start.sh $MYDIR/s
cd $MYDIR
python worker1.py $COH $SIGMA 1 
./start.sh 1 $COH $SIGMA 1 5000 &
cd ..

