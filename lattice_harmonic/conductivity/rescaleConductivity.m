function bs_moments_rescaled = rescaleConductivity(out,bs_moments)
% Rescale conductivities (& resistivities) using sum rules

%% Constants & definitions
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

        % Create parameter arrays
        pFit(gg).P    = [t*temp(gg,1) gamma(gg,1) trap(gg,1)];
        pH(gg).P      = [t*Txy(gg,1) gamma(gg,1) trap(gg,1)];
    end
end

%% Calculate sum rules/rescale factors

for hh = 1:length(out)
    tic
    % lineshapeFit    = @(f) conductivity_eval2(f,pFit(hh).P,DEPTH_ER);
    % lineshapeH      = @(f) conductivity_eval2(f,pFit(hh).P,DEPTH_ER);
    % 
    % SFit(hh)        = integral(lineshapeFit,0,inf);
    % SH(hh)          = integral(lineshapeH,0,inf);

    lineshapeFit    = conductivity_eval2(FREQ_HZ,pFit(hh).P,DEPTH_ER);
    lineshapeH      = conductivity_eval2(FREQ_HZ,pH(hh).P,DEPTH_ER);
    
    SFit(hh)        = trapz(real(lineshapeFit));
    SH(hh)          = trapz(real(lineshapeH));

    rescaleFactor(hh) = SH(hh)/SFit(hh);
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

end










