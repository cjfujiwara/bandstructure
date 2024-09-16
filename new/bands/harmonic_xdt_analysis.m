%% Caclulate band strcuture, tunneling, and wannier
% opts = struct;
% opts.doPlot = 0;
% npt=constants;
% npt.numStates=51;
% npt.numK=1001;
% npt.K=linspace(-1,1,npt.numK)';
% dK = mode(diff(npt.K));
% npt.K = [-1:dK:1]';
% if npt.K(end)~=1
%     npt.K=[npt.K; 1];
% end
% npt.numK = length(npt.K);
% npt.depth=[.1:.1:10]; 
% npt = getBandStructure(npt,opts);   % calculate band structure
% npt = calculateTunneling(npt);      % calculate tunneling elements
% showTunnelingDepth(npt)

%% Caclulate band strcuture, tunneling, and wannier
opts = struct;
opts.doPlot = 0;
npt=constants;
npt.numStates=51;
npt.numK=1001;
npt.K=linspace(-1,1,npt.numK)';
dK = mode(diff(npt.K));
npt.K = [-1:dK:1]';
if npt.K(end)~=1
    npt.K=[npt.K; 1];
end
npt.numK = length(npt.K);
npt.depth=[60:5:150]; 
npt = getBandStructure(npt,opts);   % calculate band structure
npt = calculateTunneling(npt);      % calculate tunneling elements
showTunnelingDepth(npt)

%% Caclulate band strcuture, tunneling, and wannier
opts = struct;
opts.doPlot = 1;

npt=constants;
npt.depth=[2.5]; 
npt = getBandStructure(npt,opts);   % calculate band structure
npt = calculateTunneling(npt);      % calculate tunneling elements

wannier_opts.bands = [1:4];
npt = wannier(npt,wannier_opts);           % calculate wannier function
npt = calculateWannierMoments(npt);
npt = calculateWannierHarmonicCoupling(npt);
%% Wannier Plotting
doShowWannier = 1;
doAnimateWannier = 0;

if doShowWannier
[hF_wannier] = showWannier(npt,wannier_opts);           % calculate wannier function    
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
end

%% Calculate 1D spectrum with Harmonic Confinement

% physical constants
amu = 1.66053907e-27;
h =6.626e-34;
m= 40*amu;
aL = 527e-9;

%

f_x_latt = 40;
f_y_latt = 40;
f_z_latt = 30;

fx_xdt = 42.5;
fy_xdt = 30.9;
fz_xdt = 224;

f_radial = sqrt((f_x_latt*f_y_latt) +(fx_xdt*fy_xdt)+ (f_z_latt)^2);
f_vertical = sqrt(f_x_latt^2+f_y_latt^2+fz_xdt^2);


% calculation parameters
harmonic_opts = struct;
harmonic_opts.NumSites =501;
harmonic_opts.MaxTunnelingOrder = 51;
harmonic_opts.NumBands = 4;

%42.5 ± 0.3 xdt trap freq 04/2024
%30.9 ± 0.6*



% XY Lattice
harmonic_opts.omega = 2*pi*65;
% harmonic_opts.omega=2*pi*1;
% harmonic_opts.omega = npt.fr/10;
harmonic_opts.Omega = 0.5*m*harmonic_opts.omega^2*aL^2/h;
[npt,harmonic_output_H] = calculateLatticeHarmonicSpectrum(npt,harmonic_opts);

% Z Direction
harmonic_opts.omega = 2*pi*230; % XDT Vertical trap frequency
harmonic_opts.Omega = 0.5*m*harmonic_opts.omega^2*aL^2/h;
[npt,harmonic_output_V] = calculateLatticeHarmonicSpectrum(npt,harmonic_opts);

%% Fit Lowest Band Energy to linear DOS (ie. best harmonic approximation)
harmonic_output_H = fitHOtoFirstBand(harmonic_output_H);
harmonic_output_V = fitHOtoFirstBand(harmonic_output_V);

%%

showLatticeHarmonic(harmonic_output_H,npt);
showLatticeHarmonic(harmonic_output_V,npt);

%%
showBandProjections(harmonic_output_H)
showBandProjections(harmonic_output_V)

%% 
% opts=struct;
% opts.Indeces = 'auto';
% opts.Indeces = [1 2 46 47 110 111 112 113 354 355];
% showLatticeHarmonicWavefunction(npt,harmonic_output_H,opts);

%% Thermodynamical Analysis

calculateThermodynamics(npt,harmonic_output_H,harmonic_output_H,harmonic_output_V);



 