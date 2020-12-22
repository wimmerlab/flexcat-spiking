# Simulate three single trials, for different magnitude of stimulus fluctuations (sigma).
# Save spike times rather than population firing rates.

python decision_network_sparse.py 1 8 1 spikes
mv decision_network_sparse_1.mat single_trial_sigma_8.mat
python decision_network_sparse.py 1 20 1 spikes
mv decision_network_sparse_1.mat single_trial_sigma_20.mat
python decision_network_sparse.py 1 36 1 spikes
mv decision_network_sparse_1.mat single_trial_sigma_36.mat
