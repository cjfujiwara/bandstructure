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
hf.Position(2:4) = [100 700 500];

ax=axes;
hold on
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',16,'fontname','Calibri Light');
co=get(gca,'colororder');
xlabel('magnetic field B (G)');
ylabel('energy (kHz)');
hold on
%
plot([180 230],[0 0],'-','linewidth',1,'color',[.5 .5 .5]);

m0=plot(eng_2_B_m0(E1),(E1-0*2.5)*f*1e-3,'-','color',co(2,:),'linewidth',3);
plot(eng_2_B_m0(E2),(E2-0*2.5)*f*1e-3,'-','color',co(2,:),'linewidth',3)
plot(eng_2_B_m0(E3),(E3-0*2.5)*f*1e-3,'-','color',co(2,:),'linewidth',3)
plot(eng_2_B_m0(E4),(E4-0*2.5)*f*1e-3,'-','color',co(2,:),'linewidth',3)

m1=plot(eng_2_B_m1(E1),(E1-0*2.5)*f*1e-3,':','color',co(1,:),'linewidth',3);
plot(eng_2_B_m1(E2),(E2-0*2.5)*f*1e-3,':','color',co(1,:),'linewidth',3)
plot(eng_2_B_m1(E3),(E3-0*2.5)*f*1e-3,':','color',co(1,:),'linewidth',3)
plot(eng_2_B_m1(E4),(E4-0*2.5)*f*1e-3,':','color',co(1,:),'linewidth',3)


xlim([195 208]);
ylim([-150 850]);

legend([m0 m1],{'z','x,y'},'fontsize',16,...
    'location','southeast');



% 


