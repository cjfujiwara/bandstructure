
%% Run the bootstrap on the moments 
% Analyze the 1st, 2nd, and 3rd moments. 
% Fit the 1st moment to sinuisoidal oscillation
% Fit the 2nd moment to a linear increase (for fitting)

[bs_moments,hF]=bootstrap_com(composite_data);

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
               
    for jj=1:length([bs_moments(nn).OscBootStat])
        bs_data=bs_moments(nn).OscBootStat{jj};
        C_site_bs = bs_data(:,1)*1e-6/aL;
        S_site_bs = bs_data(:,2)*1e-6/aL;
        sigma_bs = -1i*(omega./force_invsec(jj)).*(1i*C_site_bs-S_site_bs);
        rho_bs = 1./sigma_bs;
    end
    


end

%% Bootstrap frequency dependent C-S into  sigma(omega) and rho(omega)

%% Run the Bootstrap on the Spectrum
% Only do this if you really mean to, since it will take your computer a
% few hours to run
doRunBootstrap = true;
if doRunBootstrap
    for nn=1:length(bs_moments)
        out(nn)=conductivity_fit_bootstrap(bs_moments(nn));
    end
end


