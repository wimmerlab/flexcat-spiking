# flexcat-spiking - Flexible categorization dynamics in spiking neural networks of perceptual decision making

Spiking neural network simulations from the following paper:

> Prat-Ortega, G., Wimmer, K., Roxin, A. and de la Rocha, J. (2020). Flexible categorization in perceptual decision making. Nature Communications. In Press.

A complete description of the models is provided in the methods of the paper and the supplementary material.

The repository is organized into the following folders:
1. simulation: Code to run the neural network simulations. The network is implemented in [Brian2](https://briansimulator.org/) using Python 3. Installation instructions can be found [here](https://brian2.readthedocs.io/en/stable/introduction/install.html). 
2. analysis: MATLAB scripts to collect the simulation data and to compute psychometric curves and psychophysical kernels.
3. results: Simulation results (psychometric curves and psychophysical kernels) and figures.

Python scritps that reproduce all figures shown in the paper, including Fig. 5, can be found here: bitbucket.org/delaRochaLab/flexible-categorization.

If you run into any issues when running the code please contact Klaus Wimmer (wimmer.klaus@gmail.com).

## Implemented networks

We have implemented two spiking neural networks of perceptual decision making.

1. Network of leaky integrate-and fire neurons with sparse, random connectivity (probability of connection between neurons is 0.1) and current-based synapses with exponential decay. These simulations can be found in simulation/sparse_network/ and results are shown in Fig. 5 of Prat-Ortega et al., 2021.
2. Network of leaky integrate-and fire neurons with all-to-all connectivity and conductance-based synapses with AMPA, GABA and NMDA receptor dynamics. These simulations can be found in simulation/Wang2002_network/ and results are shown in Suppl. Fig. 4 of Prat-Ortega et al., 2021.

## Running simulations

The Python scripts run a single trial and save the results to a file (see start_single_trials.sh). A single simulation run will typically take a few seconds. Multiple trials, as necessary for estimating psychometric curves and psychophysical kernels, can be simulated by running the Python scripts in parallel (e.g. on a cluster). The most basic way of parallelization is to start multiple proceses, as illustrated in the included shell scripts. Note that these scripts are included for illustrative purposes and have to be adapted to the available computer infrastructure. 
