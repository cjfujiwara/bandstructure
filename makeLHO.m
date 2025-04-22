function lho = makeLHO(npt)
%MAKELHO Summary of this function goes here
%   Detailed explanation goes here

if nargin==0
    npt=struct;
end
if ~isfield(npt,'Depth_Er')
    npt.Depth_Er = 2.5;
%     npt.Depth_Er = 5;

end
if ~isfield(npt,'Nsites')
    npt.Nsites = 601;
end
if ~isfield(npt,'TrapFrequency_Hz')
    npt.TrapFrequency_Hz = 66.8;
end
if ~isfield(npt,'TunnelOrder')
    npt.TunnelOrder = 11;
end
if ~isfield(npt,'Bands')
   npt.Bands = 1;
end
%% Calculate Lattice

lattice = constants;
wannier_opts=struct;
wannier_opts.Bands = 1;

lattice.depth = npt.Depth_Er;
lattice = calculateBandStructure(lattice);   % Band structure
% lattice = calculateBandGaps(lattice);      % Band Gap
lattice = calculateTunneling(lattice);       % Tunneling
lattice = wannier(lattice,wannier_opts);     % Wannier functions
lattice = calculateWannierMoments(lattice);  % Wannier Matrix Elements

%% Caluclate Lattice + Harmonic Properties

% XY Lattice
lho_opts = struct;
lho_opts.omega              = 2*pi*npt.TrapFrequency_Hz;
lho_opts.NumSites           = npt.Nsites;
lho_opts.MaxTunnelingOrder  = npt.TunnelOrder;
lho_opts.HarmonicBands      = npt.Bands;
lho_opts.Omega              = 0.5*lattice.m*lho_opts.omega^2*...
    (lattice.lambda/2)^2/lattice.h;
[lattice,lho] = ...
    calculateLHOSpectrum_sband(lattice,lho_opts);


%% Fit to linear constant density of states
% Fit lowest band to linear dispersion
lho = fitHOtoFirstBand(lho);

%% Calculate Thermodynamics

showLHO_Eigenstates(lattice,lho)
end

