function [dB,B] = lattice_shift_79(Ux,Uy,Uz)
cc = constants;

a_bg = 166.978;
Delta = 6.910;
B0 = 202.15;

mu=cc.m/2;
a0=cc.a0;
h=cc.h;
hbar=cc.hbar;
fr=cc.fr;

R=-65*a0;
% R=-98*a0;

% Zero Point energy 1/2*hbar*omega_x+1/2*hbar*omega_y+1/2*hbar*omega_z
E0=0.5*(sqrt(4*Ux)+sqrt(4*Uy)+sqrt(4*Uz))*h*fr;




dB=Delta./(1-2*hbar^2./(mu*R*E0*a_bg*a0));
B = B0 + dB;

end

