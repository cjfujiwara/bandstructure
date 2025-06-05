%% Analayzse Density and Hubbard U
U = zeros(length(bs_moments),1);
peak_charge=zeros(length(bs_moments),1);

for bb = 1:length(composite_data)
    U(bb)=composite_data(bb).Hubbard.U(1);
    
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

%% vsU figures

hF_vs_U=figure(4);
hF_vs_U.Color='w';
clf
hF_vs_U.Name='Gamma';
hF_vs_U.Position=[50 50 400 250];
set(gcf,'color','w');

ax1=axes;
co=get(gca,'colororder');
errorbar((U/t).^2,gamma(:,1)/t/(2*pi),gamma(:,2)/t/(2*pi),...
    'o','color','k','markerfacecolor',[.5 .5 .5],...
    'linewidth',1)

xlabel('interaction strength U^2/t^2')
ylabel('current dissapation \Gamma/t')
set(gca,'fontsize',8)
xlim([0 36])
ylim([0 .05]);
hold on
%% Temperature and Density
hF_vs_U_density=figure(5);
hF_vs_U_density.Color='w';
clf
hF_vs_U_density.Name='TempDensity';
hF_vs_U_density.Position=[50 50 400 250];
set(gcf,'color','w');


ax_sub=axes;
yyaxis left
errorbar((U/t).^2,temp(:,1),temp(:,2),...
    'o','color',.5*co(1,:),'markerfacecolor',co(1,:),'markersize',6)
ylabel('spectral T/t')
set(gca,'fontsize',8)
xlim([0 36])
ylim([0 3])

yyaxis right


errorbar((U/t).^2,peak_charge(:,1)*kappa,peak_charge(:,2)*kappa,...
    'o','color',.5*co(2,:),'markerfacecolor',co(2,:))
ylabel('peak charge density n_0');
ylim([0 .2])
linkaxes([ax1 ax_sub],'x');
xlabel('interaction strength U^2/t^2')

%% Gamma over nt
hF_vsU_Gamman=figure(6);
hF_vsU_Gamman.Color='w';
hF_vsU_Gamman.Name= 'Gamma_n0';
hF_vsU_Gamman.Position=[50 50 400 250];
clf
errorbar((U/t).^2,gamma(:,1)/t/(2*pi)./peak_charge(:,1),gamma(:,2)/t/(2*pi),...
    'o','color','k','markerfacecolor',[.5 .5 .5],...
    'linewidth',1)
xlabel('U^2/t^2');
ylabel('\Gamma/(n_0 t)');
%% Temperature Comparison

hF_vsU_Temp_Compare=figure(7);
hF_vsU_Temp_Compare.Color='w';
hF_vsU_Temp_Compare.Position=[50 50 250 250];
hF_vsU_Temp_Compare.Name='TemperatureCompare';
clf
set(gcf,'color','w');

plot([0 5],[0 5],'k--');
hold on
errorbar(temp(:,1),Tx(:,1),Tx(:,2),Tx(:,2),temp(:,2),temp(:,2),...
    'o','color',co(1,:),'markerfacecolor',co(1,:));
axis equal tight
hold on
xlim([1 3]);
ylim([1 3]);
xlabel('spectral T/t');
ylabel('gaussian size T/t');

%% Central Resistivity

hF_vs_U_resistivity=figure(8);
hF_vs_U_resistivity.Color='w';
hF_vs_U_resistivity.Position=[100 100 400 250];
clf
axes
p0=errorbar((U/t).^2,rho0(:,1),rho0(:,2),...
    'o','color',.5*co(3,:),'markerfacecolor',co(3,:));
hold on
xlabel('U^2/t^2')
ylabel('\rho [1/\sigma_0]')
ylim([0 .07])