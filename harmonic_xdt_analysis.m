%% Introduction
% This script calculates the properties of an optical lattice with an
% additional harmonic confinement.

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));

a = fileparts(curpath);
addpath(a);addpath(genpath(a));

%% Initialize
% Define parameters of calculation
wannier_opts            = struct;
wannier_opts.Bands      = [1];

lattice                 = constants;
lattice.depth           = [3.5]; 
lattice.WannierBands    = wannier_opts.Bands;

%% Flags
doShowBandStructure = false;
doShowWannier = false;
doAnimateWannier = false;

%% Caclulate Band Properties
% Calculate the band structure

lattice = calculateBandStructure(lattice);   % calculate band structure

if doShowBandStructure
    hF_band = showBandStructure(lattice,wannier_opts);
end

%% Calculate Tunneling Propertiess
% Calculate the tunneling matrix elements
lattice = calculateTunneling(lattice);      % calculate tunneling elements

%% Calculate Wannier
% Calculate the wannier functions, specify which bands you want to
% calculate

lattice = wannier(lattice,wannier_opts);                % Calculate wannier function
lattice = calculateWannierMoments(lattice);             % Dipole matrix elements in wannier basis

if doShowWannier
% Show the Wannier function
    hF_wannier = showWannier(lattice,wannier_opts);           % calculate wannier function 
end
%% Make Suboutput
band_inds = [1 2 3 4 5];
xL = [-30 30];
kL = [-20 20];

i1 = find(lattice.K_extended>=kL(1),1);
i2 = find(lattice.K_extended>=kL(2),1);

i3 = find(lattice.X_extended>=kL(1),1);
i4 = find(lattice.X_extended>=kL(2),1);

maxTunnelSite = 15;

hubbard                 = struct;
hubbard.h               = lattice.h;
hubbard.hbar            = lattice.hbar;
hubbard.a0              = lattice.a0;
hubbard.lambda          = lattice.lambda;
hubbard.Er              = lattice.Er;
hubbard.fr              = lattice.fr;
hubbard.numStates       = lattice.numStates;
hubbard.numK            = lattice.numK;
hubbard.K               = lattice.K;
hubbard.depth           = lattice.depth;
hubbard.bandEigenValue  = lattice.bandEigenValue(band_inds,:,:);
hubbard.Tunneling       = lattice.Tunneling(band_inds,1:maxTunnelSite,:);

BG_1D=zeros(size(hubbard.bandEigenValue,3),1);
BG_2D=zeros(size(hubbard.bandEigenValue,3),1);
BG_3D=zeros(size(hubbard.bandEigenValue,3),1);
for kk=1:size(hubbard.bandEigenValue,3)
    Es_max = max(hubbard.bandEigenValue(1,:,kk));
    Es_min = min(hubbard.bandEigenValue(1,:,kk));
    Ep_min = min(hubbard.bandEigenValue(2,:,kk));
    BG_1D(kk) = Ep_min-Es_max;
    BG_2D(kk) = (Es_min+Ep_min)-(Es_max+Es_max);
    BG_3D(kk) = (Es_min+Es_min+Ep_min)-(Es_max+Es_max+Es_max);
end

hubbard.BandGap1D = BG_1D;
hubbard.BandGap2D = BG_2D;
hubbard.BandGap3D = BG_3D;

% lattice_out.numX            = lattice.numX;
% lattice_out.K_extended      = lattice.K_extended(i1:i2);
% lattice_out.Wannier_K       = lattice.Wannier_K(i1:i2,:,:);
% lattice_out.X_extended      = lattice.X_extended(i3:i4);
% lattice_out.Wannier_X       = lattice.Wannier_X(i3:i4,:,:);
% lattice_out.Wannier_X_Harmonic = lattice.Wannier_X_Harmonic(i3:i4,:,:);
% lattice_out.Harmonic_Length = lattice.Harmonic_Length;


%% Wannier Animation
% Animate the wannier functions if you specified different lattice depths
if doAnimateWannier
    tempfile = fullfile(tempdir,'animate.gif');
    for kk=1:length(hF_wannier)    
        frame = getframe(hF_wannier(kk));
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);  
        if kk == 1
            imwrite(A,map,tempfile,'gif','LoopCount',Inf,'DelayTime',1);
        else
            if kk==length(hF_wannier)
                imwrite(A,map,tempfile,'gif','WriteMode','append','DelayTime',1);
            else
                imwrite(A,map,tempfile,'gif','WriteMode','append','DelayTime',.1);
            end
        end        
    end
    copyfile(tempfile,'wannier.gif','f');
end

%% Wannier Harmonic Coupling
% Calculate the matrix coupling element induced from a harmonic potential
% on the wannier states (this is primarily important for multi-band
% physics)
%
% <w_m(x_i)|x^2|w_n(x_j)>
% <w_m(x_i)|x^1|w_n(x_j)>

% lattice = calculateWannierHarmonicCoupling2(lattice);
% lattice = calculateWannierHarmonicCoupling(lattice);

%% Calculate 1D spectrum with Harmonic Confinement

% Numerical Settings
Nsites = 601;
TunnelOrder = 11;
HarmonicBands = 1;

% XY Lattice
horz_opts = struct;
horz_opts.omega = 2*pi*65;
horz_opts.NumSites = Nsites;
horz_opts.MaxTunnelingOrder = TunnelOrder;
horz_opts.HarmonicBands = HarmonicBands;
horz_opts.Omega = 0.5*lattice.m*horz_opts.omega^2*(lattice.lambda/2)^2/lattice.h;
[lattice,harmonic_output_H] = calculateLHOSpectrum_sband(lattice,horz_opts);

% Z Direction
vert_opts = struct;
vert_opts.omega = 2*pi*300;
vert_opts.NumSites = Nsites;
vert_opts.MaxTunnelingOrder = TunnelOrder;
vert_opts.HarmonicBands = HarmonicBands;
vert_opts.Omega = 0.5*lattice.m*vert_opts.omega^2*(lattice.lambda/2)^2/lattice.h;
% [lattice,harmonic_output_V] = calculateLHOSpectrum(lattice,harmonic_opts);
[lattice,harmonic_output_V] = calculateLHOSpectrum_sband(lattice,vert_opts);

% Fit lowest band to linear dispersion
harmonic_output_H = fitHOtoFirstBand(harmonic_output_H);
harmonic_output_V = fitHOtoFirstBand(harmonic_output_V);

%% Show Spectrum

hF_x=showLatticeHarmonic(harmonic_output_H,lattice);
hF_z=showLatticeHarmonic(harmonic_output_V,lattice);
hF_z.Position(1) = hF_x.Position(1)+hF_x.Position(3)+5;

%% Eigen States
showLHO_Eigenstates(lattice,harmonic_output_H)
showLHO_Eigenstates(lattice,harmonic_output_V)

%% Project Eigenstates onto band original bands
% Mostly useful for multi band stuff

% showBandProjections(harmonic_output_H)
% showBandProjections(harmonic_output_V)

%% Show Eigenvectors
% Show eigenvectors and also convolve with the wannier functions
opts=struct;
opts.Indeces = 'auto';
opts.Indeces = [1 2 50 51];
% opts.Indeces = [1:100];

% showLatticeHarmonicWavefunction(lattice,harmonic_output_H,opts);

%% Thermodynamical Analysis

% calculateThermodynamics(lattice,harmonic_output_H,harmonic_output_H,harmonic_output_V);



 