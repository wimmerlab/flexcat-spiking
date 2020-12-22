function show_acc_vs_sigma(filename)
%SHOW_ACC_VS_SIGMA  Plot psychometric curve (accuracy vs. sigma vs. time)
%   SHOW_ACC_VS_SIGMA(filename) plots the results contained in the
%   specified file


load(filename,'R1_end','R2_end','acc','sigma','t','PI');
my_path = fileparts(filename);

figure
hold on
h = [];
for i=1:size(t,1)
    h(i) = plot(sigma, acc(i,:),'s-','linewidth',2);
    %plot(sigma, PC_corr(i,:),'s--','linewidth',2)
end
xlabel('Stimulus fluctuations \sigma_S (pA)')
ylabel('Accuracy')
hold on
% plot(sigma,PI,'k--','linewidth',2)
l = legend(h,string(t));
title(l,'T (s)')

figsave(gcf,fullfile(my_path, 'acc_vs_sigma_and_T'),{'fig','png'})
