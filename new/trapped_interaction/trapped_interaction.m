function [funcs] = trapped_interaction
%TRAPPED_INTERACTION Summary of this function goes here
%   Detailed explanation goes here
% Comes from the Busch paper 
%
% "Two Cold Atoms in a Harmonic Trap" Busch (1997)
%
% Assumes contact interaction of 4*pi*a_0*delta(r)
%
% Units of energy of sqrt(hbar/(m*omega)) where m is the mass of SINGLE
% particle and omega is the single particle trap frequency. See euation (1)
% and equation (2)
%
% The problem separates in a center of mass and relative motion
% hamiltonian.  The center-off-mass is easy to solve and is just the
% harmonic oscillator and so we only solve the relative motion problem
% expressed in equation (3).
%
% The energy E of the relativer motional relates to the scattering length 
% a_0 via eq. (16).
a = @(E) (sqrt(2)*gamma(-E/2+3/4)./gamma(-E/2+1/4)).^(-1);
%
% In the code we subtract off the zero point energy 1.5*hbar*omega.


Ebounds = [-20 0.5:2:10.5]; % Energy zones to consider
Np=5e5;                     % number of points
E0 = 1.5;                   % zero point energy

aMat = zeros(length(Ebounds)-1,Np);
funcs={};
for jj=1:length(Ebounds)-1
    Evec=linspace(Ebounds(jj)+1E-5,Ebounds(jj+1)-1E-5,Np);
    a_out = a(Evec);
    foo = @(a_in) interp1(a_out,Evec-E0,a_in); % solve but subtract of E0
    funcs{jj}=foo;    
    aMat(jj,:)=a_out;
end

end

