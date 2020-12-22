function [t,PK,PK_PI,n_trials] = calc_PK(my_dir, my_filename, i_max, sigma, time)
%CALC_PK  Compute psychophysical kernels
%   CALC_PK(my_dir, my_filename, i_max, sigma, time) computes the
%   psychophysical kernels for various sigmas (strength of stimulus
%   fluctuations) for different stimulus durations given by time.


smooting_window = 20;      % in ms

n_trials = zeros(1,numel(sigma));
PK = {};
PK_PI = {};
t = {};
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
    
    % decision (probability correct) read out from the neural firing rates at different times during the stimulus interval
    for i_time = 1:size(time,1)
        % decision read out from the corresponding time interval
        R1_end = mean(R1(idx,time(i_time,1):time(i_time,2)),2);
        R2_end = mean(R2(idx,time(i_time,1):time(i_time,2)),2);
        choice = (R1_end - R2_end) > 0;        % trial-to-trial choice based on the neural firing rates
        
        % decision consistent with actually shown stimulus (evidence for left
        % or right in the current trial) rather than the target valus
        choice_PI = mean(I1(:,1:time(i_time,2)),2) - mean(I2(:,1:time(i_time,2)),2) > 0;     % trial-to-trial choice for a perfect integrator
        
        
        % psychophysical kernel
        
        % prepare stimulus
        t{i_time} = 1:time(i_time,2);
        stim = I1(:,t{i_time}) - I2(:,t{i_time});
        stim = convn(stim,1/smooting_window*ones(1,smooting_window),'same');
        
        % compute AROC
        for i=1:size(stim,2)
            PK{i_time}(i_sigma,i) = calc_roc(choice(idx), stim(idx,i), 0);
            PK_PI{i_time}(i_sigma,i) = calc_roc(choice_PI(idx), stim(idx,i), 0);
        end
        
        
    end
    
end


