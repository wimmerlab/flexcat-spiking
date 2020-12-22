function collect_data(my_dir, my_filename, i_max, sigma, with_remove_files)
%COLLECT_DATA  Combine single trial simulation data in one Matlab file
%   COLLECT_DATA(my_dir, my_filename, i_max, sigma) collects the data up to
%   a maximum number of trials given by i_max and for a range of sigmas.

if nargin < 5
    with_remove_files = false;
end

for i_sigma = 1:numel(sigma)
    
    all_files = dir(sprintf('%s_%d/%s*.mat', my_dir, sigma(i_sigma), my_filename));
    fprintf('%d files found\n',numel(all_files));
    
    for i=1:min(numel(all_files),i_max)
        
        myfile = [all_files(i).folder, filesep, all_files(i).name];
        load (myfile);
        if with_remove_files
            system(sprintf('rm %s',myfile));
        end
        
        if i == 1
            R1 = NaN(i_max,numel(rateE1));          % we will always have exactly i_max trials, but they may be filled up with NaNs
            R2 = NaN(i_max,numel(rateE1));
            I1 = NaN(i_max,numel(rateE1));
            I2 = NaN(i_max,numel(rateE1));
        end
        
        I1(i,:) = I_E1;
        I2(i,:) = I_E2;
        R1(i,:) = rateE1;
        R2(i,:) = rateE2;
    end
    
    filename = [all_files(i).folder, filesep, 'sim_result.mat'];
    if exist(filename,'file')
        disp('File exists already, save anyway.')
    end
    save(filename, 'I1','I2','R1','R2','t','-v7');
    
end
