function [rhoAvg, rhoAvgRS] = averageResistivity(bs_moments,out,bs_moments_rescaled,rescaledOut)
%% Constants & definitions
t           = 563.4;    % Hz
FREQ_HZ     = 40:80;    % Hz
DEPTH_ER    = 2.5;      % recoil

freqThresh  = 0.30;     %

%% Get data
for aa = 1:length(out)
        
        % Spectral fit temperature
        pdT = fitdist(out(aa).bootstat(:,1)/t,'normal');
        temp(aa,1)=pdT.mu; 
        temp(aa,2)=pdT.sigma;
        
        % Spectral fit gamma
        pdG = fitdist(out(aa).bootstat(:,2),'normal');
        gamma(aa,1)=pdG.mu; 
        gamma(aa,2)=pdG.sigma;

        % Spectral fit trap frequency
        pdf = fitdist(out(aa).bootstat(:,3),'normal');
        trap(aa,1)=pdf.mu; 
        trap(aa,2)=pdf.sigma;

        % Create parameter arrays
        P = [t*temp(aa,1) gamma(aa,1) trap(aa,1)];

        lineshape = conductivity_eval2(FREQ_HZ,P,DEPTH_ER);

        I = find(diff(sign(imag(lineshape))));                                      % index of zero crossings
        m = (imag(lineshape(I+1))-imag(lineshape(I)))./(FREQ_HZ(I+1)-FREQ_HZ(I));   % slope
        fStar(aa,1) = -imag(lineshape(I))./m+FREQ_HZ(I);                                   
end      

%% Average resistivity
for bb = 1:length(out)
    
    % Find frequencies within X% of fStar
    iFind = find(abs(1-bs_moments(bb).Frequency_Hz/fStar(bb)) < freqThresh);

    w = 1./bs_moments(bb).rhoErr(iFind);
    rhoAvg(bb,1) = sum(real(bs_moments(bb).rho(iFind)).*w)./sum(w);
    rhoAvg(bb,1) = mean(real(bs_moments(bb).rho(iFind)));
    rhoAvg(bb,2) = sqrt(sum(real(bs_moments(bb).rhoErr(iFind)).^2))./(2*length(iFind));
end

if nargin > 2
    %% Get rescaled data
    for cc = 1:length(rescaledOut)
            
            % Spectral fit temperature
            pdT = fitdist(rescaledOut(cc).bootstat(:,1)/t,'normal');
            temp(cc,1)=pdT.mu; 
            temp(cc,2)=pdT.sigma;
            
            % Spectral fit gamma
            pdG = fitdist(rescaledOut(cc).bootstat(:,2),'normal');
            gamma(cc,1)=pdG.mu; 
            gamma(cc,2)=pdG.sigma;
    
            % Spectral fit trap frequency
            pdf = fitdist(rescaledOut(cc).bootstat(:,3),'normal');
            trap(cc,1)=pdf.mu; 
            trap(cc,2)=pdf.sigma;
    
            % Create parameter arrays
            P = [t*temp(cc,1) gamma(cc,1) trap(cc,1)];
    
            lineshape = conductivity_eval2(FREQ_HZ,P,DEPTH_ER);
    
            I = find(diff(sign(imag(lineshape))));                                      % index of zero crossings
            m = (imag(lineshape(I+1))-imag(lineshape(I)))./(FREQ_HZ(I+1)-FREQ_HZ(I));   % slope
            fStarRS(cc,1) = -imag(lineshape(I))./m+FREQ_HZ(I);                                   
    end  
    
    %% Average rescaled resistivity
    for bb = 1:length(rescaledOut)
        
        % Find frequencies within X% of fStar
        iFind = find(abs(1-bs_moments_rescaled(bb).Frequency_Hz/fStar(bb)) < freqThresh);
        
        w = 1./bs_moments_rescaled(bb).rhoErr(iFind);
        rhoAvgRS(bb,1) = sum(real(bs_moments_rescaled(bb).rho(iFind)).*w)./sum(w);

        rhoAvgRS(bb,1) = mean(real(bs_moments_rescaled(bb).rho(iFind)));
        rhoAvgRS(bb,2) = sqrt(sum(real(bs_moments_rescaled(bb).rhoErr(iFind)).^2))./(2*length(iFind));
    end
end
end