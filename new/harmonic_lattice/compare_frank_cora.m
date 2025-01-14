frank=load('frank_55Hz_2.5Er.mat');
frank=frank.frank;

harmonic_output_H;

e_fc = frank.EigenValues-frank.EigenValues(1);
e_cf = harmonic_output_H.EigenValues-harmonic_output_H.EigenValues(1);

hF=figure(30);
hF.Color='w';
clf


%
ax1=subplot(231);
co=get(gca,'colororder');
plot(e_fc,'ko','markersize',6,...
    'markerfacecolor',co(1,:));
hold on
xlabel('eigenindex');
ylabel('energy (Hz)');
grid on
title('fc energy 2.5 Er, 55 Hz');

%
ax2=subplot(234);
co=get(gca,'colororder');
plot(e_cf,...
    'ko','markersize',6,...
    'markerfacecolor',co(2,:));
hold on
xlabel('eigenindex');
ylabel('energy (Hz)');
grid on
title('cf energy 2.5 Er, 55 Hz');

linkaxes([ax1 ax2],'xy');
xlim([0 75]);
ylim([0 3000]);

%
axd1=subplot(232);
imagesc(abs(frank.Dipole));
set(gca,'ydir','normal');
hold on
title('|fc dipole|');
axis equal tight
colorbar
xlabel('eigenindex');
ylabel('eigenindex');
%
axd2=subplot(235);
imagesc(abs(harmonic_output_H.DipoleOperator));
set(gca,'ydir','normal');
hold on
title('|cf dipole|');
axis equal tight
colorbar

linkaxes([axd1 axd2],'xy','color');
xlim([1 70]);
ylim([1 70]);
caxis(axd1,[0 20]);
caxis(axd2,[0 20]);
xlabel('eigenindex');
ylabel('eigenindex');

% energy difference
ax_deng=subplot(233);
plot(e_fc(1:200)-e_cf(1:200),'ko','markersize',6,...
    'markerfacecolor',co(3,:));
hold on
xlabel('eigenindex');
ylabel('\epsilon_{fc}-\epsilon_{cf} (Hz)');
grid on
title('energy difference');
xlim([0 75]);




ax_ddipole=subplot(236);
cf_dip_1=diag(harmonic_output_H.DipoleOperator,1);
fc_dip_1=diag(frank.Dipole,1);

plot(abs(fc_dip_1),'-',...
    'color',co(1,:));
hold on
plot(abs(cf_dip_1),'-',...
    'color',co(2,:));
xlim([0 75])
xlabel('eigenindex');
ylabel('$|\langle n+1|\hat{x}|n\rangle|~(a_L)$ ','interpreter','latex',...
    'fontsize',16);
legend({'fc','cf'});
title('dipole operator');


