%% Analayzse Density and Hubbard U
% U = zeros(length(bs_moments),1);
peak_charge=zeros(length(bs_moments),1);

for bb = 1:length(composite_data)
%     U(bb)=composite_data(bb).Hubbard.U(1);
    
    n=bs_moments(bb).Density_PeakGaussCharge;
    n_me=[];
    for jj=1:length(n)
        n_me=[n_me; n{jj}];
    end

    pd_charge = fitdist(n_me,'normal');
    peak_charge(bb,1)=pd_charge.mu; 
    peak_charge(bb,2)=pd_charge.sigma;
    
end

%% Density Scalars
kappa = 0.5*(2^(-3/2)); % 0.5 for spin, 2^(3/2) for gauss
kappa = 1;


T = Tx;     % harmonic  
% T = temp;   % spectral

%% Gamma

hF_vs_T=figure(104);
hF_vs_T.Color='w';
clf
hF_vs_T.Name='Gamma';
hF_vs_T.Position=[50 50 400 250];
set(gcf,'color','w');

ax1=axes;
co=get(gca,'colororder');
errorbar(T(:,1),gamma(:,1)/t/(2*pi),gamma(:,2)/t/(2*pi),...
    'o','color','k','markerfacecolor',[.5 .5 .5],...
    'linewidth',1)

xlabel('T/t')
ylabel('current dissapation \Gamma/t')
set(gca,'fontsize',8)
xlim([0 3.5])
ylim([0 .05]);
hold on
%% Temperature and Density
hF_vs_T_density=figure(105);
hF_vs_T_density.Color='w';
clf
hF_vs_T_density.Name='TempDensity';
hF_vs_T_density.Position=[50 50 400 250];
set(gcf,'color','w');


ax_sub=axes;
yyaxis left
errorbar(T(:,1),temp(:,1),temp(:,2),...
    'o','color',.5*co(1,:),'markerfacecolor',co(1,:),'markersize',6)
ylabel('spectral T/t')
set(gca,'fontsize',8)
xlim([0 3.5])
ylim([0 3.5])

yyaxis right


errorbar(T(:,1),peak_charge(:,1)*kappa,peak_charge(:,2)*kappa,...
    'o','color',.5*co(2,:),'markerfacecolor',co(2,:))
ylabel('peak charge density n_0');
ylim([0 .2])
linkaxes([ax1 ax_sub],'x');
xlabel('T/t')

%% Gamma over nt
hF_vsT_Gamman=figure(106);
hF_vsT_Gamman.Color='w';
hF_vsT_Gamman.Name= 'Gamma_n0';
hF_vsT_Gamman.Position=[50 50 400 250];
clf
errorbar(T(:,1),gamma(:,1)/t/(2*pi)./peak_charge(:,1),gamma(:,2)/t/(2*pi),...
    'o','color','k','markerfacecolor',[.5 .5 .5],...
    'linewidth',1)
xlabel('T/t');
ylabel('\Gamma/(n_0 t)');
xlim([0 3.5])
ylim([0 .35]);
%% Temperature Comparison

hF_vsT_Temp_Compare=figure(107);
hF_vsT_Temp_Compare.Color='w';
hF_vsT_Temp_Compare.Position=[50 50 250 250];
hF_vsT_Temp_Compare.Name='TemperatureCompare';
clf
set(gcf,'color','w');

plot([0 5],[0 5],'k--');
hold on
errorbar(temp(:,1),Tx(:,1),Tx(:,2),Tx(:,2),temp(:,2),temp(:,2),...
    'o','color',co(1,:),'markerfacecolor',co(1,:));
axis equal tight
hold on
xlim([.7 3.5]);
ylim([0.7 3.5]);
xlabel('spectral T/t');
ylabel('gaussian size T/t');

%% Central Resistivity

hF_vs_T_resistivity=figure(108);
hF_vs_T_resistivity.Color='w';
hF_vs_T_resistivity.Position=[100 100 400 250];
clf
axes
p0=errorbar(T(:,1),rho0(:,1),rho0(:,2),rho0(:,2),T(:,2),T(:,2),...
    'o','color',.5*co(3,:),'markerfacecolor',co(3,:));
hold on
xlabel('T/t')
ylabel('\rho [1/\sigma_0]')
ylim([0 .12])

