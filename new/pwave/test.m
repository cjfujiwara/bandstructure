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
E1 = linspace(-1,2.49,1e5);
E2 = linspace(2.51,4.47,1e5);

figure
clf
co=get(gca,'colororder');
plot(eng_2_B_m1(E1),E1-2.5,'-','color',co(1,:),'linewidth',2)
hold on
plot(eng_2_B_m1(E2),E2-2.5,'-','color',co(1,:),'linewidth',2)
% 
plot(eng_2_B_m0(E1),E1-2.5,'--','color',co(1,:),'linewidth',2)
hold on
plot(eng_2_B_m0(E2),E2-2.5,'--','color',co(1,:),'linewidth',2)

xlim([190 210]);