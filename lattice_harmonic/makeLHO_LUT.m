function [LUT] = makeLHO_LUT(Depth_Er,TrapFrequency_Hz)

Nsites = 301;
TunnelOrder = 11;
HarmonicBands=[1];

if nargin == 0
    Depth_Er = 2.5;
    TrapFrequency_Hz = 55:.01:75;
end

if nargin ==1
    TrapFrequency_Hz = 55:.01:75;
end
%% Calculate Lattice

lattice = constants;
wannier_opts=struct;
wannier_opts.Bands = HarmonicBands;

lattice.depth = Depth_Er;
lattice = calculateBandStructure(lattice);   % Band structure
lattice = calculateTunneling(lattice);       % Tunneling
lattice = wannier(lattice,wannier_opts);     % Wannier functions
lattice = calculateWannierMoments(lattice);  % Wannier Matrix Elements

%% Make LUT


D2 = zeros(Nsites,Nsites,length(TrapFrequency_Hz));
ENG = zeros(Nsites,length(TrapFrequency_Hz));

t1=now;
for nn=1:length(TrapFrequency_Hz)
    fprintf([num2str(nn) '/' num2str(length(TrapFrequency_Hz)) ' ' num2str(TrapFrequency_Hz(nn)) ' Hz '])

    lho_opts = struct;
    lho_opts.omega              = 2*pi*TrapFrequency_Hz(nn);
    lho_opts.NumSites           = Nsites;
    lho_opts.MaxTunnelingOrder  = TunnelOrder;
    lho_opts.HarmonicBands      = HarmonicBands;   
    lho_opts.Omega              = 0.5*lattice.m*lho_opts.omega^2*...
    (lattice.lambda/2)^2/lattice.h;
    
    [~,lho] = calculateLHOSpectrum_sband(lattice,lho_opts);
    D2(:,:,nn) = lho.DipoleOperator;
    ENG(:,nn) = lho.EigenValues(:);
    disp('done');
end
t2=now;

tsec=round((t2-t1)*60*60*24,1);
disp(['LUT evaluated in ' num2str(tsec) ' sec.'])

%%

% [EE1,EE2]=meshgrid(ENG,ENG);

LUT= struct;
LUT.LATTICE                = lattice;
LUT.D2_aL                  = D2;
LUT.ENG_Hz                 = ENG;
LUT.TRAP_Hz                 = TrapFrequency_Hz;


end

