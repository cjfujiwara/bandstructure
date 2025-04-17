function output = makeLHO(npt)
%MAKELHO Summary of this function goes here
%   Detailed explanation goes here

if nargin==0
    npt=struct;
end
if ~isfield(npt,'depth_Er')
    npt.depth = 2.5;
end
if ~isfield(npt,'Nsites')
    npt.Nsites = 601;
end
if ~isfield(npt,'freq_radial_Hz')
    npt.freq_radial = 67;
end
if ~isfield(npt,'freq_axial_Hz')
    npt.freq_axial = 300;
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

lattice.depth = npt.depth;
lattice = calculateBandStructure(lattice);   % Band structure
% lattice = calculateBandGaps(lattice);      % Band Gap
lattice = calculateTunneling(lattice);       % Tunneling
lattice = wannier(lattice,wannier_opts);     % Wannier functions
lattice = calculateWannierMoments(lattice);  % Wannier Matrix Elements


%% Caluclate Lattice + Harmonic Properties

% XY Lattice
horz_opts = struct;
horz_opts.omega = 2*pi*npt.freq_radial;
horz_opts.NumSites = npt.Nsites;
horz_opts.MaxTunnelingOrder = npt.TunnelOrder;
horz_opts.HarmonicBands = npt.Bands;
horz_opts.Omega = 0.5*lattice.m*horz_opts.omega^2*(lattice.lambda/2)^2/lattice.h;
[lattice,harmonic_output_H] = calculateLHOSpectrum_sband(lattice,horz_opts);

% Z Direction
vert_opts = struct;
vert_opts.omega = 2*pi*npt.freq_axial;
vert_opts.NumSites = npt.Nsites;
vert_opts.MaxTunnelingOrder = npt.TunnelOrder;
vert_opts.HarmonicBands = npt.Bands;
vert_opts.Omega = 0.5*lattice.m*vert_opts.omega^2*(lattice.lambda/2)^2/lattice.h;
% [lattice,harmonic_output_V] = calculateLHOSpectrum(lattice,harmonic_opts);
[lattice,harmonic_output_V] = calculateLHOSpectrum_sband(lattice,vert_opts);


%% Fit to linear constant density of states
% Fit lowest band to linear dispersion
harmonic_output_H = fitHOtoFirstBand(harmonic_output_H);
harmonic_output_V = fitHOtoFirstBand(harmonic_output_V);

%% Calculate Thermodynamics
% calculateThermodynamics(lattice,harmonic_output_H,harmonic_output_H,harmonic_output_V);


showLHO_Eigenstates(lattice,harmonic_output_H)
showLHO_Eigenstates(lattice,harmonic_output_V)
end

