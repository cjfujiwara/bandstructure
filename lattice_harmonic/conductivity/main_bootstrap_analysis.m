
%% Run the bootstrap on the moments 
% Analyze the 1st, 2nd, and 3rd moments. 
% Fit the 1st moment to sinuisoidal oscillation
% Fit the 2nd moment to a linear increase (for fitting)

doRunBootstrap_moments = true;
if doRunBootstrap_moments
    [bs_moments,hF]=bootstrap_com(composite_data);
end

%% Get Conductivity
omega_xdt   = 2*pi*42;          % [1/s] XDT Trap Frequency
amu         = 1.66054e-27;      % [kg] atomic mass unit
m           = 40*amu;           % [kg] potassium-40 mass
aL          = 532e-9;           % [m] lattice spacing
umperv      = 3.6;              % [um/V] piezo conversion factor
hbar        = 1.05457182e-34;   % [Js] reduce planck's constant

for nn=1:length(bs_moments)
    v2=zeros(length(bs_moments(nn).Params),1);
    for jj=1:length([bs_moments(nn).Params])
        v2(jj) = unique([bs_moments(nn).Params{jj}.conductivity_ODT2_mod_amp]);
    end
    omega           = 2*pi*[bs_moments(nn).Frequency_Hz];omega=omega(:);
    x0              = (v2*umperv*1e-6);
    force_invsec    = (m*omega_xdt^2*x0*aL)/hbar;
    bs_moments(nn).force_invsec ...
                    = force_invsec;
    C_site          =  [bs_moments(nn).C_um]*1e-6/aL;C_site=C_site(:);
    S_site          =  [bs_moments(nn).S_um]*1e-6/aL;S_site=S_site(:);
    CErr_site       =  abs([bs_moments(nn).CErr_um]*1e-6/aL);CErr_site=CErr_site(:);
    SErr_site       =  abs([bs_moments(nn).SErr_um]*1e-6/aL);SErr_site=SErr_site(:);
    bs_moments(nn).sigma ...
                    = -1i*(omega./force_invsec).*(-1i*C_site+S_site);
    bs_moments(nn).sigmaErr ...
                    = 1i*(omega./force_invsec).*(-1i*CErr_site+SErr_site);
    bs_moments(nn).rho ...
                    = (-1i*(omega./force_invsec).*(-1i*C_site+S_site)).^-1;
    bs_moments(nn).rhoErr ...
                    = -(bs_moments(nn).sigmaErr./abs(bs_moments(nn).sigma.^2));                  
end

%% Average real resistivity data near omega_star or real conductivity data peak
doAverageResistivity = true;
if doAverageResistivity
    if (exist('bs_moments_Gibbs_rescaled','var') && exist('bs_moments_rescaled','var'))
        [rhoAvg, rhoAvgRS, rhoAvgGRS] = averageResistivity(bs_moments,bs_moments_rescaled,bs_moments_Gibbs_rescaled);  
    elseif exist('bs_moments_rescaled','var')
        [rhoAvg, rhoAvgRS] = averageResistivity(bs_moments,bs_moments_rescaled);
    else
        [rhoAvg] = averageResistivity(bs_moments);
    end
end
%% Run the Bootstrap on the Spectrum
% Only do this if you really mean to, since it will take your computer a
% few hours to run
doRunBootstrap = true;
if doRunBootstrap
    for nn=1:length(bs_moments)
        out(nn)=conductivity_fit_bootstrap(bs_moments(nn));
    end
end

%% Rescale conductivities
doRescale = true;
if doRescale
    [bs_moments_rescaled,bs_moments_Gibbs_rescaled,TG] = rescaleConductivity(out,bs_moments,composite_data);
end

%% Run the bootstrap on the rescaled spectrum
% Only do this if you really mean to, since it will take your computer a
% few hours to run
doRunRescaledBootstrap = true;
if doRunRescaledBootstrap
    for nn=1:length(bs_moments_rescaled)
        rescaledOut(nn)=conductivity_fit_bootstrap(bs_moments_rescaled(nn));
    end
end
%% Run the bootstrap on the Gibbs rescaled spectrum
% Only do this if you really mean to, since it will take your computer a
% few hours to run
doRunGibbsRescaledBootstrap = false;
if doRunGibbsRescaledBootstrap
    for nn=1:length(bs_moments_rescaled)
        GibbsRescaledOut(nn)=conductivity_fit_bootstrap(bs_moments_Gibbs_rescaled(nn));
    end
end

