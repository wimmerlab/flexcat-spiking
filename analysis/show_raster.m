function show_raster(file, max_t)
%SHOW_RASTER  Show single trial example.
%   SHOW_RASTER(file) plots the raster plot and stimulus input for a single
%   trial.

load(file);

spks_t_E1 = spks_t_E1(:)';
spks_i_E1 = spks_i_E1(:)';
spks_t_E2 = spks_t_E2(:)';
spks_i_E2 = spks_i_E2(:)';

figure('position',[ 440    94   560   704]);

subplot(11,1,1:5);
plot(reshape([spks_t_E1; spks_t_E1; NaN(size(spks_t_E1))],1,[]),...
    reshape( max(spks_i_E1) + [spks_i_E1-0.5; spks_i_E1+0.5; NaN(size(spks_t_E1))],1,[]), ...
    'color','r','linestyle','-','linewidth',2);
hold on
plot(reshape([spks_t_E2; spks_t_E2; NaN(size(spks_t_E2))],1,[]),...
    reshape( 0 + [spks_i_E2-0.5; spks_i_E2+0.5; NaN(size(spks_t_E2))],1,[]), ...
    'color','b','linestyle','-','linewidth',2);
axis([0 1000*max_t 0 2*max(spks_i_E1)]);
axis off

subplot(11,1,6:8);
plot(t/1000,rateE2,'b')
hold on
plot(t/1000,rateE1,'r')
axis([0 max_t 0 100]);
ylabel('Rate (Hz)')
set(gca,'xtick',0:2:max_t)
set(gca,'xticklabel',[])
box off
set(gca,'fontsize',12)

subplot(11,1,9:11);
plot(t/1000,I_E2 * 1000,'b','linewidth',1)
hold on
plot(t/1000,I_E1 * 1000,'r','linewidth',1)
axis([0 max_t -30 60]);
xlabel('Time (s)')
set(gca,'xtick',0:2:max_t)
ylabel({'Stimulus','strength (pA)'})
box off
set(gca,'fontsize',12)
set(gca,'ytick',([0 50]));
