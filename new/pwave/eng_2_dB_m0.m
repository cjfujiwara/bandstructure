function [B] = eng_2_dB_m0(V,eng)
a0 = 5.29177210903E-11;
amu = 1.66053906660E-27;
h = 6.62607004E-34;
hbar = h/(2*pi);

m = 39.96399848 * amu;
mu = m/2;

Delta = -19.89;
B0 = 198.803;
Vbg = -(108*a0)^3;
R = 49.4*a0;

dB = Delta * (1 - (Vbg*(1./V-2*mu*eng/(R*hbar^2))).^(-1)).^(-1);

B = dB + B0;
end

