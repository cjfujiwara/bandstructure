
% Fundamental constants
a0 = 5.29177210903E-11;
amu = 1.66053906660E-27;
h = 6.62607004E-34;
hbar = h/(2*pi);
m = 39.96399848 * amu;
mu = m/2;

E1 = linspace(-3,2.48,1e5);
E2 = linspace(2.501,4.475,1e5);
E3 = linspace(4.5,6.4,1e5);
E4 = linspace(6.5,8.8,1e5);

%%
loc = 'G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Composite P-wave';

f100 =  'data_100Er.mat';
f200 = 'data_200Er_new.mat';
f200 =  'data_200Er.mat';

f300 =  'data_300Er.mat';
%  f300 = 'data_300Er_old.mat';
% f300 = 'data_265Er.mat';

d1 = load(fullfile(loc,f100));
d2 = load(fullfile(loc,f200));
d3 = load(fullfile(loc,f300));
% 
% d2New = load(fullfile(loc,f200New));
% d2b= load(fullfile(loc,f200));


hf1=figure(201);
clf
hf1.Color='w';
hf1.Position(2:4) = [-200 1000 600];
co=get(gca,'colororder');

fs = [85.188e3 122.531E3 151.545e3];
% fs = [85.188e3 122.531E3 140e3];

clear ps
clear ps2

theory1=1;
theory0=1;
unitarity=1;

for ii=1:length(fs)  
    % Trap freqency in Hz
    f = fs(ii);
    
    % useful parameters
    omega = 2*pi*f;             % Trap frequency in rad
    E0 = hbar*omega;            % Trap energy
    d = sqrt(hbar/(mu*omega));  % Oscillator lenght

    % Energy [normalized] to volume [m^3]
    eng2V_m1 = @(eng) d^3./eq34(eng,1);
    eng2V_m0 = @(eng) d^3./eq33(eng,1);

    % Energy [normalized] to field
    eng_2_B_m1 = @(eng) eng_2_dB_m1(eng2V_m1(eng),eng*hbar*omega);
    eng_2_B_m0 = @(eng) eng_2_dB_m0(eng2V_m0(eng),eng*hbar*omega);

    if theory1
    % m1
    plot(eng_2_B_m1(E1),(E1-2.5)*f*1e-3,'-','color',co(ii,:),'linewidth',2);    
    hold on
    ps(ii)=plot(eng_2_B_m1(E2),(E2-2.5)*f*1e-3,'-','color',co(ii,:),'linewidth',2);
    end
    if theory0
    % m0
    ps2(ii)=plot(eng_2_B_m0(E1),(E1-2.5)*f*1e-3,'--','color',co(ii,:),'linewidth',2);    
    hold on
    plot(eng_2_B_m0(E2),(E2-2.5)*f*1e-3,'--','color',co(ii,:),'linewidth',2);
    end
    if  unitarity
        y = ps(ii).YData;
        j=find(y>f*1e-3,1);
        x = ps(ii).XData(j);
        plot(x,f*1e-3*1,'s','color',.5*co(ii,:),'linewidth',2,...
            'markerfacecolor',co(ii,:),'markersize',15);

    end
end

hold on

p1=errorbar(d1.data_process.B,2*d1.data_process.df,...
    d1.data_process.df_err,d1.data_process.df_err,...
    d1.data_process.B_err,d1.data_process.B_err,...
    'o','markerfacecolor',co(1,:),...
    'markeredgecolor',.5*co(1,:),'color',co(1,:)*.5,...
    'linewidth',2,'markersize',8); 
hold on

p2=errorbar(d2.data_process.B,2*d2.data_process.df,...
    d2.data_process.df_err,d2.data_process.df_err,...
    d2.data_process.B_err,d2.data_process.B_err,...
    'o','markerfacecolor',co(2,:),...
    'markeredgecolor',.5*co(2,:),'color',co(2,:)*.5,...
    'linewidth',2,'markersize',8); 
hold on

p3=errorbar(d3.data_process.B,2*d3.data_process.df,...
    d3.data_process.df_err,d3.data_process.df_err,...
    d3.data_process.B_err,d3.data_process.B_err,...
    'o','markerfacecolor',co(3,:),...
    'markeredgecolor',.5*co(3,:),'color',co(3,:)*.5,...
    'linewidth',2,'markersize',8); 
hold on

xlabel('magnetic field (G)');
ylabel('energy (kHz)');

ylim([-80 155]);
xlim([199.2 201.5]);
% xlim([round(min(all_B)-0.15,1) round(max(all_B)+0.15,1)]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',16,'xminorgrid','on','yminorgrid','on','fontname','calibri light');

% s1 = '$100~E_\mathrm{R}~(\Delta_\mathrm{O}=85.19~\mathrm{kHz})$';
% s2 = '$200~E_\mathrm{R}~(\Delta_\mathrm{O}=122.54~\mathrm{kHz})$';
% s3 = '$300~E_\mathrm{R}~(\Delta_\mathrm{O}=151.15~\mathrm{kHz})$';
s1 = '$\Delta_\mathrm{O}=85.19~\mathrm{kHz}$';
s2 = '$\Delta_\mathrm{O}=122.54~\mathrm{kHz}$';
s3 = '$\Delta_\mathrm{O}=151.15~\mathrm{kHz}$';
legend([p1 p2 p3],{s1,s2,s3},'location','northwest',...
    'interpreter','latex','fontsize',10);

if theory1
legend([p1 p2 p3 ps ],{s1,s2,s3, 'x,y','x,y','x,y'},'location','northwest',...
    'interpreter','latex','fontsize',10);
end

if theory0
legend([p1 p2 p3 ps ps2],{s1,s2,s3, 'x,y','x,y','x,y','z','z','z'},'location','northwest',...
    'interpreter','latex','fontsize',10);
end