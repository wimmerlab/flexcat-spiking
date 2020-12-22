function show_PK(filename)
%SHOW_PK  Plot psychophysical kernels.
%   SHOW_PK(filename) plots the results contained in the
%   specified file

load(filename,'sigma_PK','time','t','PK');
my_path = fileparts(filename);


% show results
my_sigma = sigma_PK;
my_time = time(end:-1:1);

figure('position',[192   449   944   259]);

for i_my_sigma = 1:numel(my_sigma)
    i_sigma = find(my_sigma(i_my_sigma) == sigma_PK);
    subplot(1,3,i_my_sigma); hold on
    axis([0 max(time/1000) 0.45 0.75])
    plot(xlim,[0.5 0.5],'k--')
    for i_my_time = 1:numel(my_time)
        i_time = find(my_time(i_my_time) == time);
        plot(t{i_time}/1000,PK{i_time}(i_sigma,:),'linewidth',2);
    end
    title(sprintf('sigma = %d',my_sigma(i_my_sigma)))
end



figsave(gcf,fullfile(my_path, 'PK_vs_sigma_and_T'),{'fig','png'})

