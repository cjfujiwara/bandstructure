function out=buchler_interaction
%BUCHLER_INTERACTION Summary of this function goes here
%   Detailed explanation goes here
%PhysRevLett.104.090402
%
% Equation (8) gives the onsite interaction term for the hubbard model
% single band
%
%
% U = W*w0^2/(1+W*gamma)
%
% W = 8/pi*Er*a_s/a_L
%
% Where a_L is the lattice spacing and Er=hbar^2*k_L^2/(2*m) where kL is
% the lattice light wavevector.

% In terms of the harmonic oscillator length scale this becomes ...
%
% W = hbar*omega_HO*0.5*u0^(-3/4)*a_s/a_HO

str = 'PhysRevLett.104.090402';

w02 = [2.412 5.954 9.483 12.63 15.50];  % |w_0|^2;  unitless
gamma = [0.6 2.3 3.5 4 2];              % gamma; in units of 1/Er
u0 = [4 8 12 16 20];                    % in units of Er


% gamma in units not of 1/Er but in 1/hbar*omega
gamma_HO = gamma.*sqrt(4*u0);

gamma_HO = gamma_HO;

W_HO = @(u0,as) 0.5*u0.^(-3/4).*as;

W_HO = @(u0,as) 2^(7/2)/pi^2.*(4*u0).^(-3/4).*as;

W_HO = @(u0,as) 2^(7/2)/pi^2.*(4*u0).^(-3/4).*as;



out.U_04 = @(as) w02(1).*W_HO(u0(1),as)./(1+W_HO(u0(1),as)*gamma_HO(1));
out.U_08 = @(as) w02(2).*W_HO(u0(2),as)./(1+W_HO(u0(2),as)*gamma_HO(2));
out.U_12 = @(as) w02(3).*W_HO(u0(3),as)./(1+W_HO(u0(3),as)*gamma_HO(3));
out.U_16 = @(as) w02(4).*W_HO(u0(4),as)./(1+W_HO(u0(4),as)*gamma_HO(4));
out.U_20 = @(as) w02(5).*W_HO(u0(5),as)./(1+W_HO(u0(5),as)*gamma_HO(5));


end

