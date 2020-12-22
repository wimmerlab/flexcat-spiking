# Simulate three single trials, for different magnitude of stimulus fluctuations (sigma).
# Save spike times rather than population firing rates.

python decision_network_all_to_all.py 8 40 1 spikes
mv decision_network_all_to_all_1.mat single_trial_sigma_40.mat
python decision_network_all_to_all.py 8 80 1 spikes
mv decision_network_all_to_all_1.mat single_trial_sigma_80.mat
python decision_network_all_to_all.py 8 140 1 spikes
mv decision_network_all_to_all_1.mat single_trial_sigma_140.mat
