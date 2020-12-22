'''

Decision network as in:

Wang, X.-J. 
Probabilistic decision making by slow reverberation in cortical circuits. 
Neuron, 2002, 36, 955-968.

@author: Klaus Wimmer
wimmer.klaus@gmail.com

-----------------------------------------------------------------------------
In this variant, the stimulus input is implemented as a fluctuating current 
injected into excitatory neurons of the coding populations (in contrast to 
Wang 2002, where the stimulus input is modeled as input from an external 
population with Poisson firing rates).
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

# stimulus and simulation parameters
coh = icoh                       # coherence of random dots
sigma = isigma / 100. * 10 * pA  # standard deviation of stimulus input (stimulus fluctations)
Imu0 = 10 * pA                   # injected current to both population (stimulus input at zero coherence)
Imu1 = 10 * pA                   # selective stimulus input at highest coherence
tau_OU = 100 * ms                # time constant of stimulus fluctuations
stim_on = 1000 * ms              # stimulus onset
stim_off = 21000 * ms            # stimulus offset
runtime = 21000 * ms             # total simulation time

# external noise inputs
N_ext = 1000                     # number of external Poisson neurons  
rate_ext_E = 2400 * Hz / N_ext   # external Poisson rate for excitatory population
rate_ext_I = 2400 * Hz / N_ext   # external Poisson rate for inhibitory population

# network parameters
N = 2000                         # number of neurons
f_inh = 0.2                      # fraction of inhibitory neurons
NE = int(N * (1.0 - f_inh))      # number of excitatory neurons (1600)
NI = int(N * f_inh)              # number of inhibitory neurons (400)
fE = 0.15                        # coding fraction
subN = int(fE * NE)              # number of neurons in decision pools (240)

# neuron parameters
El = -70.0 * mV                  # resting potential
Vt = -50.0 * mV                  # firing threshold
Vr = -55.0 * mV                  # reset potential
CmE = 0.5 * nF                   # membrane capacitance for pyramidal cells (excitatory neurons)
CmI = 0.2 * nF                   # membrane capacitance for interneurons (inhibitory neurons)
gLeakE = 25.0 * nS               # membrane leak conductance of excitatory neurons
gLeakI = 20.0 * nS               # membrane leak conductance of inhibitory neurons
refE = 2.0 * ms                  # refractory periodof excitatory neurons
refI = 1.0 * ms                  # refractory period of inhibitory neurons

# synapse parameters
V_E = 0. * mV                    # reversal potential for excitatory synapses
V_I = -70. * mV                  # reversal potential for inhibitory synapses
tau_AMPA = 2.0 * ms              # AMPA synapse decay
tau_NMDA_rise = 2.0 * ms         # NMDA synapse rise
tau_NMDA_decay = 100.0 * ms      # NMDA synapse decay
tau_GABA = 5.0 * ms              # GABA synapse decay
alpha = 0.5 * kHz                # saturation of NMDA channels at high presynaptic firing rates
C = 1 * mmole                    # extracellular magnesium concentration

# synaptic conductances
gextE = 2.1 * nS                 # external -> excitatory neurons (AMPA)
gextI = 1.62 * nS                # external -> inhibitory neurons (AMPA)
gEEA = 0.05 * nS / NE * 1600     # excitatory -> excitatory neurons (AMPA)
gEIA = 0.04 * nS / NE * 1600     # excitatory -> inhibitory neurons (AMPA)
gEEN = 0.165 * nS / NE * 1600    # excitatory -> excitatory neurons (NMDA)
gEIN = 0.13 * nS / NE * 1600     # excitatory -> inhibitory neurons (NMDA)
gIE = 1.3 * nS / NI * 400        # inhibitory -> excitatory neurons (GABA)
gII = 1.0 * nS / NI * 400        # inhibitory -> inhibitory neurons (GABA)

# synaptic footprints
Jp = 1.7                         # relative synaptic strength inside a selective population (1.0: no potentiation))
Jm = 1.0 - fE * (Jp - 1.0) / (1.0 - fE)

# neuron equations
# note the "(unless refractory)" statement serves to clamp the membrane voltage during the refractory period;
# otherwise the membrane potential continues to be integrated but no spikes are emitted
eqsE = """
   dV/dt = (- gLeakE * (V - El) - I_AMPA - I_NMDA - I_GABA - I_AMPA_ext + I_input) / CmE : volt (unless refractory)

   I_AMPA = s_AMPA * (V - V_E) : amp
   ds_AMPA / dt = - s_AMPA / tau_AMPA : siemens

   I_NMDA = s_NMDA_tot * (V - V_E) / ( 1 + exp(-0.062 * V/mvolt) * (C/mmole / 3.57) ) : amp
   s_NMDA_tot : siemens

   I_GABA = s_GABA * (V - V_I) : amp
   ds_GABA / dt = - s_GABA / tau_GABA : siemens

   I_AMPA_ext = s_AMPA_ext * (V - V_E) : amp
   ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : siemens

   I_input = stimulus(t, int(i>=subN) + int(i>=2*subN)): amp

   ds_NMDA / dt = - s_NMDA / tau_NMDA_decay + alpha * x * (1 - s_NMDA) : 1
   dx / dt = - x / tau_NMDA_rise : 1
"""

eqsI = """
   dV/dt = (- gLeakI * (V - El) - I_AMPA - I_NMDA - I_GABA - I_AMPA_ext) / CmI : volt (unless refractory)

   I_AMPA = s_AMPA * (V - V_E) : amp
   ds_AMPA / dt = - s_AMPA / tau_AMPA : siemens

   I_NMDA = s_NMDA_tot * (V - V_E) / ( 1 + exp(-0.062 * V/mvolt) * (C/mmole / 3.57) ): amp
   s_NMDA_tot : siemens

   I_GABA = s_GABA * (V - V_I) : amp
   ds_GABA / dt = - s_GABA / tau_GABA : siemens

   I_AMPA_ext = s_AMPA_ext * (V - V_E) : amp
   ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : siemens
"""

# neuron populations
popE = NeuronGroup(NE, model = eqsE, threshold = 'V > Vt', reset = 'V = Vr', refractory = refE, method='euler')
popI = NeuronGroup(NI, model = eqsI, threshold = 'V > Vt', reset = 'V = Vr', refractory = refI, method='euler')
popE1 = popE[:subN]
popE2 = popE[subN:2*subN]
popE3 = popE[2*subN:]

# connections involving AMPA synapses

# excitatory -> excitatory connection through AMPAR
C_EE_AMPA = Synapses(popE, popE, 'w : siemens', on_pre = 's_AMPA += w', delay = 0.5 * ms, method='euler')
C_EE_AMPA.connect()
C_EE_AMPA.w[:] = gEEA
C_EE_AMPA.w[C_EE_AMPA.indices[:subN, :subN]] *= Jp
C_EE_AMPA.w[C_EE_AMPA.indices[subN:2*subN, subN:2*subN]] *= Jp
C_EE_AMPA.w[C_EE_AMPA.indices[:subN, subN:2*subN]] *= Jm
C_EE_AMPA.w[C_EE_AMPA.indices[subN:2*subN, :subN]] *= Jm
C_EE_AMPA.w[C_EE_AMPA.indices[2*subN:, :subN]] *= Jm
C_EE_AMPA.w[C_EE_AMPA.indices[2*subN:, subN:2*subN]] *= Jm

# excitatory -> inhibitory connection through AMPAR
C_EI_AMPA = Synapses(popE, popI, on_pre = 's_AMPA += gEIA', delay = 0.5 * ms, method='euler')
C_EI_AMPA.connect()   # perhaps use C_DE_DI_AMPA.connect('i != j')

# external inputs (fixed background firing rates)
extinputE = PoissonInput(popE, 's_AMPA_ext', N_ext, rate_ext_E, gextE)
extinputI = PoissonInput(popI, 's_AMPA_ext', N_ext, rate_ext_I, gextI)

# connections involving NMDA synapses
C_EE_NMDA = Synapses(popE, popE, on_pre='x_pre += 1', delay = 0.5 * ms, method='euler')
C_EE_NMDA.connect(j='i')

s_NMDA_array = asarray(popE.s_NMDA)               # remark: using the asarray notation brings a massive performance improvement
s_NMDA_totE_array = asarray(popE.s_NMDA_tot)      # remark: using asarray with the subgroup (e.g. with "decisionE2" does not work, 
s_NMDA_totI_array = asarray(popI.s_NMDA_tot)      #         so we get the right neurons later in "update_NMDA" 

# calculate NMDA contributions in each time step
# contribution of NMDA can be calculated efficiently due to the all-to-all connectivity        
@network_operation(when='start')
def update_NMDA():

     sE1 = sum(s_NMDA_array[:subN])
     sE2 = sum(s_NMDA_array[subN:2*subN])
     sE3 = sum(s_NMDA_array[2*subN:])
     s_NMDA_totE_array[:subN] = gEEN * (Jp * sE1 + Jm * sE2 + Jm * sE3)
     s_NMDA_totE_array[subN:2*subN] = gEEN * (Jm * sE1 + Jp * sE2 + Jm * sE3)
     s_NMDA_totE_array[2*subN:] = gEEN * (sE1 + sE2 + sE3)
     s_NMDA_totI_array[:] = gEIN * (sE1 + sE2 + sE3)

# connections involving GABA synapses

# inhibitory-excitatory connections
C_IE = Synapses(popI, popE, on_pre = 's_GABA += gIE', delay = 0.5 * ms, method='euler')
C_IE.connect()

# inhibitory-inhibitory connections
C_II = Synapses(popI, popI, on_pre = 's_GABA += gII', delay = 0.5 * ms, method='euler')
C_II.connect()

# set initial conditions
popE.s_NMDA_tot = gEEN * tau_NMDA_decay * 10 * Hz * 0.2
popI.s_NMDA_tot = gEIN * tau_NMDA_decay * 10 * Hz * 0.2
popE.V = Vt - 2 * mV
popI.V = Vt - 2 * mV

# stimulus input
n_pre_stim_intervals = int(np.ceil(stim_on / defaultclock.dt))
stim_duration = stim_off - stim_on
n_post_stim_intervals = int(np.ceil((runtime - stim_off) / defaultclock.dt))
I_E1 = Quantity(r_[zeros(n_pre_stim_intervals), mk_OU_stim(stim_duration = stim_duration, sigma = sigma, mu0 = Imu0 + coh / 100.0 * Imu1, tau = tau_OU, dt = defaultclock.dt), zeros(n_post_stim_intervals)],Imu0.dim)
I_E2 = Quantity(r_[zeros(n_pre_stim_intervals), mk_OU_stim(stim_duration = stim_duration, sigma = sigma, mu0 = Imu0 - coh / 100.0 * Imu1, tau = tau_OU, dt = defaultclock.dt), zeros(n_post_stim_intervals)],Imu0.dim)
I_E3 = Quantity(zeros(int(ceil(runtime / defaultclock.dt))),Imu0.dim)
stimulus = TimedArray(Quantity(np.transpose([I_E1,I_E2,I_E3]),dim = I_E1.dim), dt = defaultclock.dt)


# -----------------------------------------------------------------------------------------------
# run the simulation
# -----------------------------------------------------------------------------------------------

# record population activity
R1 = PopulationRateMonitor(popE1)
R2 = PopulationRateMonitor(popE2)

if with_save_spikes:
    # record spikes of excitatory neurons in the decision encoding populations
    SME1 = SpikeMonitor(popE1, record=True)
    SME2 = SpikeMonitor(popE2, record=True)

# run the simulation
run(runtime, report='stdout')

# smooth firing rate and downsample to 1 ms resolution
rateE1 = R1.smooth_rate(window = 'flat', width = 1 * ms) / Hz
rateE2 = R2.smooth_rate(window = 'flat', width = 1 * ms) / Hz
t = R1.t / ms
spacing = int(1 * ms / defaultclock.dt)
rateE1 = rateE1[0::spacing]
rateE2 = rateE2[0::spacing]
t = t[0::spacing]
I_E1 = I_E1 / nA
I_E2 = I_E2 / nA
I_E1 = I_E1[0::spacing]
I_E2 = I_E2[0::spacing]

if with_save_spikes:
   savemat(resfile, {'rateE1': rateE1, 'rateE2': rateE2, 't': t, 'I_E1': I_E1, 'I_E2': I_E2,
      'spks_t_E1':SME1.t / ms, 'spks_i_E1':array(SME1.i), 
      'spks_t_E2':SME2.t / ms, 'spks_i_E2':array(SME2.i)}, do_compression=True, oned_as='column')
else:
   savemat(resfile, {'rateE1': rateE1, 'rateE2': rateE2, 't': t, 'I_E1': I_E1, 'I_E2': I_E2}, 
      do_compression=True, oned_as='column')
