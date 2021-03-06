function out = constants

out=struct;

%% Physical Constants
h = 6.62607004E-34;                              % Planck's constant [J s]
hbar = h/(2*pi);                            % Reduced Planck's constant [J s]
a0=5.29177210903E-11;                       % Bohr Radius

out.h=h;
out.hbar=hbar;
out.a0=a0;
%% Atomic Parameters
m = 39.96399817 * 1.66053907E-27;           % 40K mass [kg]
out.m=m;
%% Lattice Parameters
lambda = 1053.6E-9;                         % Lattice light wavelength [m]
k_lambda =(2*pi)/lambda;                    % Lattice light wavevector
Er = hbar^2*k_lambda^2/(2*m);               % Recoil Energy [J]
fr = Er/h;                                  % Recoil Frequency [Hz]

out.lambda=lambda;
out.k_lambda=k_lambda;
out.Er=Er;
out.fr=fr;

%% Numerical parameters
out.numStates=201;
out.numK=51;

% out.numStates=25;
% out.numK=51;

out.K=linspace(-1,1,out.numK)';

end

