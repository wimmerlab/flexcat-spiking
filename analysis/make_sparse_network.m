% ------------------------------------------------------------------------
% directories
% ------------------------------------------------------------------------

sim_dir = '../simulation/sparse_network/';        % directory containing simulations
res_dir = '../results/sparse_network/';           % directory where results will be saved


% ------------------------------------------------------------------------
% (1) run simulations in Python
%
% use shell scripts to run simulations on a cluster
% e.g. start_all_coh_0.sh
%
% for single trial simulations use start_single_trials.sh
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% (2) collect the results
% ------------------------------------------------------------------------

n_trials = 5000;            % number of trials for psychometric curve (accuracy vs. sigma)
n_trials_PK = 5000;         % number of trials for PKs
sigma = 0:4:52;             % simulated sigmas
sigma_PK = [8 18 36];       % simulated sigmas for PKs

collect_data([sim_dir 'coh_1_sigma'],'worker',n_trials,sigma);
collect_data([sim_dir 'coh_0_sigma'],'worker',n_trials_PK,sigma_PK);


% ------------------------------------------------------------------------
% (3) figure showing performance vs. sigma
% ------------------------------------------------------------------------

% spike count window for computing the decision
pre_stim_time = 500;
time(1,:) = pre_stim_time + [ 900 1000 ];
time(2,:) = pre_stim_time + [ 1900 2000 ];
time(3,:) = pre_stim_time + [ 2900 3000 ];
time(4,:) = pre_stim_time + [ 3900 4000 ];
time(5,:) = pre_stim_time + [ 4900 5000 ];
time(6,:) = pre_stim_time + [ 5900 6000 ];

[R1_end,R2_end,acc,sigma,PI] = calc_acc_vs_sigma([sim_dir 'coh_1_sigma'],'sim_result',n_trials,sigma,time);

% save results
t = (time(:,end)-pre_stim_time) / 1000;
sigma = sigma / 100 * 25;       % sigma in pA; note that for simulations with coh = 1 (mu = 0.015) the average selective input to each population is 0.3750 pA
save([res_dir 'acc_vs_sigma_and_T.mat'],'R1_end','R2_end','acc','sigma','t','PI','-v7');

% plot accuracy vs. sigma
show_acc_vs_sigma([res_dir 'acc_vs_sigma_and_T.mat']);


% ------------------------------------------------------------------------
% (4) figure showing PKs for different sigmas
% ------------------------------------------------------------------------

% spike count window for computing the decision
pre_stim_time = 500;
time(1,:) = pre_stim_time + [ 1900 2000 ];
time(2,:) = pre_stim_time + [ 3900 4000 ];
time(3,:) = pre_stim_time + [ 5900 6000 ];

[t,PK,PK_PI,n_trials] = calc_PK([sim_dir 'coh_0_sigma'],'sim_result',n_trials_PK,sigma_PK,time);

% save results
time = (time(:,end)-pre_stim_time);
sigma_PK = sigma_PK / 100 * 25;       % sigma in pA; note that for simulations with coh = 1 (mu = 0.015) the average selective input to each population is 0.3750 pA
save([res_dir 'PK_vs_sigma_and_T.mat'],'sigma_PK','time','t','PK','PK_PI','n_trials_PK','-v7');

% plot PKs
show_PK([res_dir 'PK_vs_sigma_and_T.mat']);


% ------------------------------------------------------------------------
% (5) single trial raster plots
% ------------------------------------------------------------------------

show_raster([sim_dir 'single_trial_sigma_8.mat'],4.5); figsave(gcf,[res_dir 'single_trial_sigma_8'],{'fig','png'});
show_raster([sim_dir 'single_trial_sigma_20.mat'],4.5); figsave(gcf,[res_dir 'single_trial_sigma_20'],{'fig','png'});
show_raster([sim_dir 'single_trial_sigma_36.mat'],4.5); figsave(gcf,[res_dir 'single_trial_sigma_36'],{'fig','png'});
