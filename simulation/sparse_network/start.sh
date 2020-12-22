# Run many trials (indexed from I1 to I2) for a given coherence and sigma.

WORKER=$1
COH=$2
SIGMA=$3
I1=$4
I2=$5

for i in $(eval echo {$I1..$I2})
do
    ( python worker$WORKER.py $COH $SIGMA $i )
done

echo All simulations done.
