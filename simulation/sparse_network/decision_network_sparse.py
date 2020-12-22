'''

Sparse decision network from:

Prat-Ortega, G., Wimmer, K., Roxin, A. and de la Rocha, J. 
Flexible categorization in perceptual decision making.
Nature Communications, 2020, In press.

@author: Klaus Wimmer
wimmer.klaus@gmail.com

-----------------------------------------------------------------------------
This network is similar to the one in Roxin & Ledberg 2008, but with exponential 
synapses instead of instantaneous synapses.
-----------------------------------------------------------------------------

tested with:
Python 3.7.3, Brian2 2.2.2.1
Python 3.8.2, Brian2 2.3

'''

from brian2 import *
import sys
import os
from scipy.io import savemat
from scipy.signal import lfilter

# EFFICIENCY: use automatic code generation
# Note that a different random number generator is used when using C++ code generation
# Note that code generation cannot be used together with @network_operation and there is no error message for that.
(mydir,_) = os.path.splitext(__file__)  # use local directory for running cpp code
set_device('cpp_standalone',directory=mydir, debug=False)


def mk_OU_stim(stim_duration = 1000 * ms, sigma = 0.08 * nA / 5.0, mu0 = 0.40 * nA / 5.0, tau = 20.0 * ms, dt = 0.01 * ms):
# simulate UO process in discrete time => AR(1) process
# i has variance sigma^2 and mean mu

    lamda = 1.0 / (tau / dt)
    timesteps = int(stim_duration / dt)
    a = np.exp(-lamda)
    b = mu0 * (1.0 - a)
    c = sigma * np.sqrt(1 - a * a)

    i = lfilter(np.ones(1),[1, -a], b + c * np.random.randn(timesteps))

    return Quantity(i, mu0.dim)


# -----------------------------------------------------------------------------------------------
# parse command line arguments
# -----------------------------------------------------------------------------------------------

if len(sys.argv) < 4:
   print ('Need to specify: coherence, sigma and trial number')
   sys.exit("Wrong number of arguments.")

icoh = int(sys.argv[1])
isigma = int(sys.argv[2])
itrial = int(sys.argv[3])
(resfile,_) = os.path.splitext(__file__)  # results will be saved in file with the same name as the script, but with extension '.mat'
resfile += '_' + str(itrial) + '.mat'

# save spike trains in addition to population rate
with_save_spikes = len(sys.argv) == 5 and sys.argv[4] == 'spikes'

print ('coh: ', icoh, 'sigma: ', isigma, ' trial: ', itrial)


# -----------------------------------------------------------------------------------------------
# set up the simulation
# -----------------------------------------------------------------------------------------------

# stimulus parameters
coh = 1.5 * icoh                 # coherence of random dots
sigma = isigma / 100. * 25 * pA  # standard deviation of stimulus input (stimulus fluctations)
Imu0 = 25 * pA                   # injected current to both population (stimulus input at zero coherence)
Imu1 = 25 * pA                   # selective stimulus input at highest coherence
tau_OU = 20 * ms                 # time constant of stimulus fluctuations

# simulation parameters
synapse_type = 'exponential'     # choose between 'instantaneous' or 'exponential' synapses
stim_duration  = 6000 * ms       # stimulus duration
runtime = 6500 * ms              # total simulation time
random_seed_net = None           # to make simulations reproducible, specific random seed here
random_seed_run = None           # to make simulations reproducible, specific random seed here

# external noise inputs
Jext = 0.2 * mV                  # EPSP weight from external input to E and I neurons
N_ext = 1000                     # number of external Poisson neurons
rate_ext_E = 5000 * Hz / N_ext   # external Poisson rate for excitatory population
rate_ext_I = 9000 * Hz / N_ext   # external Poisson rate for inhibitory population

# network parameters
NE = 1000                        # number of neurons in each of the two excitatory populations
NI = 500                         # number of inhibitory neurons

# neuron parameters
El = -70.0 * mV                  # resting potential
Vt = -50.0 * mV                  # firing threshold
Vr = -60.0 * mV                  # reset potential
taum_e = 20.0 * ms               # membrane time constant of excitatory neurons
taum_i = 10.0 * ms               # membrane time constant of inhibitory neurons
Cm = 0.25 * nF                   # membrane capacitance

# synapse parameters
taus_e = 12.5 * ms               # time constant of excitatory exponential synapses
taus_i = 1.0 * ms                # time constant of inhibitory exponential synapses 
delay_e = 'rand() * 10 * ms'     # synaptic delays for excitatory synapses (uniform distribution with mean 5 ms)
delay_i = 'rand() * 2 * ms'      # synaptic delays for inhibitory synapses (uniform distribution with mean 1 ms)

# recurrent connectivity parameters
Jee = 0.1 * mV                   # EPSP weight from E <- E neurons (supercritical: 0.1 * mV; subcritical: 0.16 * mV)
Jie = 0.1 * mV                   # EPSP weight from I <- E neurons
Jei = -0.2 * mV                  # EPSP weight from E <- I neurons (supercritical: -0.2 * mV; subcritical: -0.22 * mV)
Cee = 100                        # number of connection from E <- E
Cie = 50                         # number of connections from I <- E
Cei = 50                         # number of connections from E <- I

# neuron equations
if synapse_type == 'instantaneous':
		  
   eqsE = """
	      dV/dt = - (V - El) / taum_e + (I / Cm) : volt
	      I = stimulus(t, int(i>=NE)): amp
	      """

   eqsI = """
	      dV/dt = - (V - El) / taum_i : volt
	      """

elif synapse_type == 'exponential':

   eqsE = """
	      dV/dt = (ge + gi - (V - El)) / taum_e + (I / Cm) : volt
	      I = stimulus(t, int(i>NE)): amp
	      dge/dt = - ge / taus_e : volt
	      dgi/dt = - gi / taus_i : volt
	      """

   eqsI = """
	      dV/dt = (ge + gi - (V - El)) / taum_i : volt
	      dge/dt = -ge / taus_e : volt
	      dgi/dt = -gi / taus_i : volt
	      """

# initialize
np.random.seed(random_seed_net)
seed(random_seed_net)            # sets the seed of Brian random number generator for stand-alone code generation

# neuron populations
decisionE = NeuronGroup(2 * NE, model = eqsE, threshold = 'V > Vt', reset='V = Vr', method = 'euler')
decisionI = NeuronGroup(NI, model = eqsI, threshold = 'V > Vt', reset = 'V = Vr', method = 'euler')
decisionE1 = decisionE[:NE]
decisionE2 = decisionE[NE:]

# external inputs

# fixed background inputs to inhibitory neurons
extinput_I = PoissonInput(decisionI, 'V', N_ext, rate_ext_I, Jext)
extinput_E = PoissonInput(decisionE, 'V', N_ext, rate_ext_E, Jext)

# recurrent connectivity

# excitatory to excitatory connections
if synapse_type == 'instantaneous':
   ConectE1E1 = Synapses(decisionE1, decisionE1, on_pre = 'V += Jee')
elif synapse_type == 'exponential':
   ConectE1E1 = Synapses(decisionE1, decisionE1, on_pre = 'ge += Jee * taum_e/taus_e')
ConectE1E1.connect(p = 1.0 * Cee / NE)
ConectE1E1.delay = delay_e

if synapse_type == 'instantaneous':
   ConectE2E2 = Synapses(decisionE2, decisionE2, on_pre = 'V += Jee')
elif synapse_type == 'exponential':
   ConectE2E2 = Synapses(decisionE2, decisionE2, on_pre = 'ge += Jee * taum_e/taus_e')
ConectE2E2.connect(p = 1.0 * Cee / NE)
ConectE2E2.delay = delay_e

# excitatory to inhibitory connections
if synapse_type == 'instantaneous':
   ConectIE1 = Synapses(decisionE, decisionI, on_pre = 'V += Jie')
elif synapse_type == 'exponential':
   ConectIE1 = Synapses(decisionE, decisionI, on_pre = 'ge += Jie * taum_i/taus_e')
ConectIE1.connect(p = 1.0 * Cie / NE)
ConectIE1.delay = delay_e

# inhibitory to excitatory connections
if synapse_type == 'instantaneous':
   ConectE1I = Synapses(decisionI, decisionE, on_pre='V += Jei')
elif synapse_type == 'exponential':
   ConectE1I = Synapses(decisionI, decisionE, on_pre='gi += Jei*taum_e/taus_i')
ConectE1I.connect(p = 1.0 * Cei / NI)
ConectE1I.delay = delay_i

# initialize random number generator
np.random.seed(random_seed_run)
seed(random_seed_run)            # sets the seed of Brian random number generator for stand-alone code generation

# set random initial conditions
decisionE.V = 'Vr + rand() * (Vt - Vr)'
decisionI.V = 'Vr + rand() * (Vt - Vr)'
if synapse_type == 'exponential':
   decisionE.ge = 0 * mV
   decisionE.gi = 0 * mV
   decisionI.ge = 0 * mV
   decisionI.gi = 0 * mV

# stimulus input
I_E1 = mk_OU_stim(stim_duration = stim_duration, sigma = sigma, mu0 = Imu0 + coh / 100.0 * Imu1, tau = tau_OU, dt = defaultclock.dt)
I_E2 = mk_OU_stim(stim_duration = stim_duration, sigma = sigma, mu0 = Imu0 - coh / 100.0 * Imu1, tau = tau_OU, dt = defaultclock.dt)

add_zeros = int(runtime / defaultclock.dt) - int(stim_duration/defaultclock.dt)
I_E1 = Quantity(np.concatenate([np.zeros(add_zeros) * nA, np.asarray(I_E1)]), dim = I_E1.dim)
I_E2 = Quantity(np.concatenate([np.zeros(add_zeros) * nA, np.asarray(I_E2)]), dim = I_E2.dim)

stimulus = TimedArray(Quantity(np.transpose([I_E1,I_E2]),dim = I_E1.dim), dt = defaultclock.dt)


# -----------------------------------------------------------------------------------------------
# run the simulation
# -----------------------------------------------------------------------------------------------

# record population activity
R1 = PopulationRateMonitor(decisionE1)
R2 = PopulationRateMonitor(decisionE2)
RI = PopulationRateMonitor(decisionI)

if with_save_spikes:
    # record spikes
   SME1 = SpikeMonitor(decisionE1, record=True) 
   SME2 = SpikeMonitor(decisionE2, record=True)
   SMI = SpikeMonitor(decisionI, record=True)

# run the simulation
run(runtime, report='stdout')

# smooth firing rate and downsample to 1 ms resolution
rateE1 = R1.smooth_rate(window = 'flat', width = 1 * ms) / Hz
rateE2 = R2.smooth_rate(window = 'flat', width = 1 * ms) / Hz
rateI = RI.smooth_rate(window = 'flat', width = 1 * ms) / Hz
t = R1.t / ms
spacing = int(1 * ms / defaultclock.dt)
rateE1 = rateE1[0::spacing]
rateE2 = rateE2[0::spacing]
rateI = rateI[0::spacing]
t = t[0::spacing]
I_E1 = I_E1 / nA
I_E2 = I_E2 / nA
I_E1 = I_E1[0::spacing]
I_E2 = I_E2[0::spacing]

if with_save_spikes:
   savemat(resfile, {'rateE1': rateE1, 'rateE2': rateE2, 'rateI': rateI, 't': t, 'I_E1': I_E1, 'I_E2': I_E2,
      'spks_t_E1':SME1.t / ms, 'spks_i_E1':array(SME1.i), 
      'spks_t_E2':SME2.t / ms, 'spks_i_E2':array(SME2.i), 
      'spks_t_I':SMI.t / ms, 'spks_i_I':array(SMI.i)}, do_compression=True, oned_as='column')
else:
   savemat(resfile, {'rateE1': rateE1, 'rateE2': rateE2, 'rateI': rateI, 't': t, 'I_E1': I_E1, 'I_E2': I_E2}, 
      do_compression=True, oned_as='column')
