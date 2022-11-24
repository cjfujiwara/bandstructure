f = 122.5E3; % Trap freqency in Hz

% Fundamental constants
a0 = 5.29177210903E-11;
amu = 1.66053906660E-27;
h = 6.62607004E-34;
hbar = h/(2*pi);
m = 39.96399848 * amu;
mu = m/2;

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

Eall = linspace(-10,10,1e5);
E1 = linspace(-3,2.4945,1e5);
E2 = linspace(2.501,4.475,1e5);
E3 = linspace(4.5,6.4,1e5);
E4 = linspace(6.5,8.8,1e5);

hf=figure(1);
clf
hf.Color='w';
hf.Position(1:4) = [10 10 300 700];
set(groot,'defaultAxesTickLabelInterpreter','latex');
ax=axes;
hold on
set(gca,'xgrid','off','ygrid','off','box','on','linewidth',3,...
    'fontsize',30,'fontname','Calibri Light','visible','on','yaxislocation','right');
co=get(gca,'colororder');

set(gca,'XTickLabel',{},'YTickLabel',{});
set(gca,'YTick',[0 5/2*f 9/2*f]*1e-3);

hold on

plot([180 235],5/2*[1 1]*f*1e-3,'k--','linewidth',2)
plot([180 235],9/2*[1 1]*f*1e-3,'k--','linewidth',2)

plot([180 230],[0 0],'-','linewidth',2,'color','k');

m1=plot(eng_2_B_m1(E1),(E1-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5);
m1.Color=m1.Color*.5;
m1=plot(eng_2_B_m1(E2),(E2-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5)
m1.Color=m1.Color*.75;
m1=plot(eng_2_B_m1(E3),(E3-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5)
m1=plot(eng_2_B_m1(E4),(E4-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5)

xlim([195 208]);
ylim([-150 620]);

hf=figure(2);
clf
hf.Color='w';
hf.Position(2:4) = [100 500 500];
set(groot,'defaultAxesTickLabelInterpreter','latex');
ax=axes;
hold on
set(gca,'xgrid','off','ygrid','off','box','on','linewidth',3,...
    'fontsize',30,'fontname','Calibri Light','visible','off','yaxislocation','right');
co=get(gca,'colororder');

set(gca,'XTickLabel',{},'YTickLabel',{});
set(gca,'YTick',[0 5/2*f 9/2*f]*1e-3);

hold on

plot([180 235],5/2*[1 1]*f*1e-3,'k--','linewidth',2)
plot([180 235],9/2*[1 1]*f*1e-3,'k--','linewidth',2)

plot([180 230],[0 0],'-','linewidth',2,'color','k');

m1=plot(eng_2_B_m1(E1),(E1-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5);
m1.Color=m1.Color*.5;
m1=plot(eng_2_B_m1(E2),(E2-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5)
m1.Color=m1.Color*.75;
m1=plot(eng_2_B_m1(E3),(E3-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5)
m1=plot(eng_2_B_m1(E4),(E4-0*2.5)*f*1e-3,'-','color','magenta','linewidth',5)

xlim([199.0567 201.0567]);
ylim([217.9167 387.9167]);

