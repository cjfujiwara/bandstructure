frank=load('frank_55Hz_2.5Er.mat');
frank=frank.frank;

harmonic_output_H;

e_fc = frank.EigenValues-frank.EigenValues(1);
dip_fc = frank.Dipole;
dip_fc2 = dip_fc/1i;
e_cf = harmonic_output_H.EigenValues-harmonic_output_H.EigenValues(1);
dip_cf = harmonic_output_H.DipoleOperator;


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

%%
imagesc(abs(dip_cf));set(gca,'ydir','normal');xlim([1 20]);ylim([1 20]);caxis([0 20])

bb=ctranspose(dip_cf).*dip_cf;
%% Lineshape Comparison

T_list = [1 2 3 4 5];
T_list = linspace(1,6,10);
t=567;
Gamma = 2*pi*15;
f_list = 10:.1:150;

% energy cf
N=200;
[exx_cf,eyy_cf]=meshgrid(e_cf(1:N),e_cf(1:N));
df_cf = eyy_cf-exx_cf;
% dipole cf
x_cf2 = ctranspose(dip_cf).*dip_cf;
x_cf2=x_cf2(1:N,1:N);
% sigma cf
sigma_cf = zeros(length(f_list),length(T_list));

% energy fc
[exx_fc,eyy_fc]=meshgrid(e_fc,e_fc);
df_fc = eyy_fc-exx_fc;
% dipole cf
x_fc2 = abs(dip_fc).^2;
% sigma fc
sigma_fc = zeros(length(f_list),length(T_list));


for tt=1:length(T_list)
    T= T_list(tt)*t;

    % CORA
    Z_cf = sum(exp(-e_cf(1:N)./T));
    fp_cf = exp(-eyy_cf/T)/Z_cf;
    fpp_cf = exp(-exx_cf/T)/Z_cf;
    A_cf = (fp_cf-fpp_cf).*x_cf2;

    % frank
    Z_fc = sum(exp(-e_fc./T));
    fp_fc = exp(-eyy_fc/T)/Z_fc;
    fpp_fc = exp(-exx_fc/T)/Z_fc;
    A_fc = (fp_fc-fpp_fc).*x_fc2;


    for ff=1:length(f_list)
        f = f_list(ff);

        f_mat_cf = 2*pi*(f - df_cf)+1i*Gamma/2;       
        C_cf = A_cf./f_mat_cf;
        sigma_cf(ff,tt) = -1i*(2*pi*f)*sum(C_cf,'all');

        f_mat_fc = 2*pi*(f - df_fc)+1i*Gamma/2;       
        C_fc = A_fc./f_mat_fc;
        sigma_fc(ff,tt) = -1i*(2*pi*f)*sum(C_fc,'all');
    end


end

figure(20)
co=jet(length(T_list));
clf

ax1=subplot(211);
ax2=subplot(212);

legStr={};
clear ps
for tt=1:size(sigma_cf,2)
    axes(ax1);
    ps(tt)=plot(f_list,real(sigma_cf(:,tt)),'-','color',co(tt,:));
    hold on
    plot(f_list,real(sigma_fc(:,tt)),'--','color',co(tt,:));
    xlabel('frequency (Hz)');
    ylabel('Re(\sigma)/\sigma_0')
title(['\Gamma=' num2str(round(Gamma,1)) ' 1/s'])

    legStr{tt}=[num2str(T_list(tt)) 't'];
    hold on
    axes(ax2);
    plot(f_list,imag(sigma_cf(:,tt)),'-','color',co(tt,:));
    hold on
    plot(f_list,imag(sigma_fc(:,tt)),'--','color',co(tt,:));
    xlabel('frequency (Hz)');
    ylabel('Im(\sigma)/\sigma_0')
title(['\Gamma=' num2str(round(Gamma,1)) ' 1/s'])

end
legend(ax1,ps,legStr);
%%
[sigma_max_cf,ind]=max(real(sigma_cf),[],1);
f_max_cf=f_list(ind);

[sigma_max_fc,ind]=max(real(sigma_fc),[],1);
f_max_fc=f_list(ind);

figure(10);
clf
 plot(T_list,f_max_cf)
 hold on
  plot(T_list,f_max_fc)


figure(21);
clf


subplot(121)
imagesc(T_list,f_list,real(sigma_cf))
xlabel('temperature (t)')
ylabel('Re(\sigma)')
ylim([30 60]);

subplot(122)
imagesc(T_list,f_list,real(sigma_fc))
xlabel('temperature (t)')
ylabel('Re(\sigma)')
ylim([30 60]);
