function [R1_end,R2_end,acc,sigma,PI] = calc_acc_vs_sigma(my_dir, my_filename, i_max, sigma, time)
%CALC_ACC_VS_SIGMA  Compute psychometric curve (accuracy vs. sigma)
%   CALC_ACC_VS_SIGMA(my_dir, my_filename, i_max, sigma, time) computes the
%   probability of correct choice vs strength of stimulus fluctuations for
%   different stimulus durations given by time.

n_trials = zeros(1,numel(sigma));
PI = zeros(1,numel(sigma));
acc = zeros(size(time,1),numel(sigma));
PC_corr = zeros(size(time,1),numel(sigma));
for i_sigma = 1:numel(sigma)
    
    all_files = dir(sprintf('%s_%d/%s*.mat', my_dir, sigma(i_sigma), my_filename));
    fprintf('%d files found\n',numel(all_files));
    
    i = 1;
    data = load ([all_files(i).folder, filesep, all_files(i).name]);
    
    maxtime = max(time(:));
    n_t = size(data.R1,2);
    
    % use zero-padding if the saved results are shorter than maxtime
    R1 = [zeros(i_max,maxtime-n_t) data.R1(1:i_max,1:min(n_t,maxtime))];
    R2 = [zeros(i_max,maxtime-n_t) data.R2(1:i_max,1:min(n_t,maxtime))];
    I1 = [zeros(i_max,maxtime-n_t) data.I1(1:i_max,1:min(n_t,maxtime))];
    I2 = [zeros(i_max,maxtime-n_t) data.I2(1:i_max,1:min(n_t,maxtime))];
    
    % idx are the valid trials
    idx = ~isnan(sum(R1,2));
    n_trials(i_sigma) = sum(idx);
    R1 = R1(idx,:);
    R2 = R2(idx,:);
    I1 = I1(idx,:);
    I2 = I2(idx,:);
    
    PI(i_sigma) = mean ( (sum(I1,2) > sum(I2,2)));     % fraction of trials where the target answer is consistent with integral of stimulus ("perfect integrator")
    
    % decision (probability correct) read out from the neural firing rates at different times during the stimulus interval
    for i_time = 1:size(time,1)
        % correct decision with respect to underlying target value (always
        % I1 > I2 in the simulations)
        R1_end{i_sigma}(i_time,:) = mean(R1(idx,time(i_time,1):time(i_time,2)),2);
        R2_end{i_sigma}(i_time,:) = mean(R2(idx,time(i_time,1):time(i_time,2)),2);
        d_r = (R1_end{i_sigma}(i_time,:) - R2_end{i_sigma}(i_time,:)) > 0;        % trial-to-trial choice based on the neural firing rates
        acc(i_time,i_sigma) = mean (d_r);
        
        % decision consistent with actually shown stimulus (evidence for left
        % or right in the current trial) rather than the target valus
        d_i = (mean(I1(idx,:),2) - mean(I2(idx,:),2)) > 0;        % trial-to-trial choice based on the stimulus
        PC_corr(i_time,i_sigma) = mean ( d_r(:) == d_i(:) );
    end
    
end

