function [bs_moments_rescaled,bs_moments_Gibbs_rescaled,TG] = rescaleConductivity(out,bs_moments,cd)
% Rescale conductivities (& resistivities) using sum rules

%% Constants & definitions
constants_conductivity;
t = 563.4; % Hz
DEPTH_ER = 2.5; % recoil

FREQ_HZ = 0:1e3;

%% Get data
% alreadyDone = true;
alreadyDone = false;

if ~alreadyDone
    clear temp
    clear gamma
    clear trap
    clear Tx
    clear Ty
    clear Txy

    for gg = 1:length(out)
        
        % Spectral fit temperature
        pdT = fitdist(out(gg).bootstat(:,1)/t,'normal');
        temp(gg,1)=pdT.mu; 
        temp(gg,2)=pdT.sigma;
        
        % Spectral fit gamma
        pdG = fitdist(out(gg).bootstat(:,2),'normal');
        gamma(gg,1)=pdG.mu; 
        gamma(gg,2)=pdG.sigma;

        % Spectral fit trap frequency
        pdf = fitdist(out(gg).bootstat(:,3),'normal');
        trap(gg,1)=pdf.mu; 
        trap(gg,2)=pdf.sigma;
        
        % Harmonic temperature X
        pdTx = fitdist(out(gg).bootstat(:,6)/t,'normal');
        Tx(gg,1)=pdTx.mu; 
        Tx(gg,2)=pdTx.sigma;

        % Harmonic temperature Y
        pdTy = fitdist(out(gg).bootstat(:,7)/t,'normal');
        Ty(gg,1)=pdTy.mu; 
        Ty(gg,2)=pdTy.sigma;

        % Calculate geometric mean harmonic temperature
        Txy(gg,1) = sqrt(Tx(gg,1).*Ty(gg,1));
        Txy(gg,2) = (Ty(gg,1).*Tx(gg,2)+Tx(gg,1).*Ty(gg,2))./(2*sqrt(Tx(gg,1).*Ty(gg,1)));

        % do Gibbs fit with bootstrap size
        Gibbs_opts = struct;
        Gibbs_opts.doGibbsRefit = 1;
        Gibbs_opts.TrapOmega = 2*pi*[trap(gg,1) trap(gg,1)];
        [cd,GibbsTemperature] = calculateGibbsTemperature(cd,Gibbs_opts);
    
        % Assign Gibbs fit properties
        TG(gg,1)                = mean(GibbsTemperature(gg).T)*kB/t*1e-9/h;
        TG(gg,2)                = std(GibbsTemperature(gg).T)*kB/t*1e-9/h;
        fig5.Gnpeak_singlon     = mean(GibbsTemperature(gg).npeak_singlon);
        fig5.Gnpeak_singlonErr  = std(GibbsTemperature(gg).npeak_singlon);
        fig5.Gnpeak_doublon     = mean(GibbsTemperature(gg).npeak_doublon);
        fig5.Gnpeak_doublonErr  = std(GibbsTemperature(gg).npeak_doublon);
        nG(gg,1)                = mean(GibbsTemperature(gg).npeak_singlon+GibbsTemperature(gg).npeak_doublon);
        nG(gg,2)                = std(GibbsTemperature(gg).npeak_singlon+GibbsTemperature(gg).npeak_doublon);

        % Create parameter arrays
        pFit(gg).P  = [t*temp(gg,1) gamma(gg,1) trap(gg,1)];
        pH(gg).P    = [t*Txy(gg,1) gamma(gg,1) trap(gg,1)];
        pG(gg).P    = [t*TG(gg,1) gamma(gg,1) trap(gg,1)];
    end
end

%% Calculate sum rules/rescale factors

for hh = 1:length(out)
    tic
    lineshapeFit    = conductivity_eval2(FREQ_HZ,pFit(hh).P,DEPTH_ER);
    lineshapeH      = conductivity_eval2(FREQ_HZ,pH(hh).P,DEPTH_ER);
    lineshapeG      = conductivity_eval2(FREQ_HZ,pG(hh).P,DEPTH_ER);
    
    SFit(hh)        = trapz(real(lineshapeFit));
    SH(hh)          = trapz(real(lineshapeH));
    SG(hh)          = trapz(real(lineshapeG));

    rescaleFactor(hh) = SH(hh)/SFit(hh);
    rescaleFactorG(hh) = SG(hh)/SFit(hh);
    toc
end

%% Rescale conductivities and resistivities

% Initialize rescaled moments structure
bs_moments_rescaled = bs_moments;

for ii = 1:length(out)
    bs_moments_rescaled(ii).sigma = bs_moments(ii).sigma.*rescaleFactor(ii);
    bs_moments_rescaled(ii).sigmaErr = bs_moments(ii).sigmaErr.*rescaleFactor(ii);
    bs_moments_rescaled(ii).rho = bs_moments(ii).rho./rescaleFactor(ii);
    bs_moments_rescaled(ii).rhoErr = bs_moments(ii).rhoErr./rescaleFactor(ii);
    bs_moments_rescaled(ii).rescaleFactor = rescaleFactor(ii);
end

%% Rescale conductivities and resistivities

% Initialize rescaled moments structure
bs_moments_Gibbs_rescaled = bs_moments;

for ii = 1:length(out)
    
    bs_moments_Gibbs_rescaled(ii).sigma = bs_moments(ii).sigma.*rescaleFactorG(ii);
    bs_moments_Gibbs_rescaled(ii).sigmaErr = bs_moments(ii).sigmaErr.*rescaleFactorG(ii);
    bs_moments_Gibbs_rescaled(ii).rho = bs_moments(ii).rho./rescaleFactorG(ii);
    bs_moments_Gibbs_rescaled(ii).rhoErr = bs_moments(ii).rhoErr./rescaleFactorG(ii);
    bs_moments_Gibbs_rescaled(ii).rescaleFactor = rescaleFactorG(ii);
end
end










