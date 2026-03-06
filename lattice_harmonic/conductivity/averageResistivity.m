function [rhoAvg,rhoAvgRS,rhoAvgGRS] = averageResistivity(bs_moments,bs_moments_rescaled,bs_moments_Gibbs_rescaled)
if nargin < 3
    rhoAvgGRS   = [];
elseif nargin < 2
    rhoAvgRS    = [];
    rhoAvgGRS   = [];
end
%% Constants & definitions
t           = 563.4;    % Hz
FREQ_HZ     = 40:80;    % Hz
DEPTH_ER    = 2.5;      % recoil

freqThresh  = 0.6;      
freqThreshL = 0.4;
freqThreshU = 1.5;

%% Get data
for aa = 1:length(bs_moments)
        
        % % Spectral fit temperature
        % pdT = fitdist(out(aa).bootstat(:,1)/t,'normal');
        % temp(aa,1)=pdT.mu; 
        % temp(aa,2)=pdT.sigma;
        % 
        % % Spectral fit gamma
        % pdG = fitdist(out(aa).bootstat(:,2),'normal');
        % gamma(aa,1)=pdG.mu; 
        % gamma(aa,2)=pdG.sigma;
        % 
        % % Spectral fit trap frequency
        % pdf = fitdist(out(aa).bootstat(:,3),'normal');
        % trap(aa,1)=pdf.mu; 
        % trap(aa,2)=pdf.sigma;
        % 
        % % Create parameter arrays
        % P = [t*temp(aa,1) gamma(aa,1) trap(aa,1)];
        % 
        % lineshape = conductivity_eval2(FREQ_HZ,P,DEPTH_ER);
        % 
        % I = find(diff(sign(imag(lineshape))));                                      % index of zero crossings
        % m = (imag(lineshape(I+1))-imag(lineshape(I)))./(FREQ_HZ(I+1)-FREQ_HZ(I));   % slope
        % fStar(aa,1) = -imag(lineshape(I))./m+FREQ_HZ(I);
        
        [M,IMAX] = max(real(bs_moments(aa).sigma));
        fStar(aa,1) = bs_moments(aa).Frequency_Hz(IMAX);
end      

%% Average resistivity
clear rhoAvg
for bb = 1:length(bs_moments)
    
    % Find frequencies within X% of fStar
    iFind = find(abs(1-bs_moments(bb).Frequency_Hz/fStar(bb)) < freqThresh);
    % iFind = find(bs_moments(bb).Frequency_Hz/fStar(bb) > freqThreshL & bs_moments(bb).Frequency_Hz/fStar(bb) < freqThreshU);
    % % iFind = find(bs_moments(bb).Frequency_Hz < 80 & bs_moments(bb).Frequency_Hz > 50);
    % iFind = find((real(bs_moments(bb).sigma).^2+abs(imag(bs_moments(bb).sigma).^2)) > 60);
    % iFind = find(sqrt((bs_moments(bb).C_um.^2+bs_moments(bb).S_um.^2)) > 0.9);

    N = length(iFind);

    % Real resistivities in frequency range
    rho = real(bs_moments(bb).rho(iFind));
    
    % Uncertainties
%     drho = abs(real(bs_moments(bb).rhoErr(iFind)/2));
    
    %Propagate uncertainty
    re = real(bs_moments(bb).sigma(iFind));
    im = imag(bs_moments(bb).sigma(iFind));
    rerr = real(bs_moments(bb).sigmaErr(iFind)/2); %divide by 2 to get 67% confidence
    ierr = imag(bs_moments(bb).sigmaErr(iFind)/2); %divide by 2 to get 67% confidence
    rhoErr = sqrt(rerr.^2.*(im.^2-re.^2).^2 + 4*re.^2.*im.^2.*ierr.^2)./(re.^2+im.^2).^2 + 1j*sqrt(ierr.^2.*(im.^2-re.^2).^2 + 4*re.^2.*im.^2.*rerr.^2)./(re.^2+im.^2).^2;
    
    % Uncertainties in Re[rho] (67% confidence)
    drho = abs(real(rhoErr));

    % Weights
    w = drho.^(-2);

    wrho  = struct;
    for loop = 1:length(rho)
        wrho(loop).rho = rho(loop);
        wrho(loop).w = w(loop);
    end
    
    %Define weighted average function
    rhoAvg_bs = bootstrp(2000,@rhoAvgfn,wrho);
    rhoAvg_pd = fitdist(rhoAvg_bs,'Normal');
    % keyboard
    %Output bootstrap weighted avg rho with bs dist sigma width as uncertainty
    rhoAvg(bb,1) = rhoAvg_pd.mu;
    rhoAvg(bb,2) = rhoAvg_pd.sigma;
    
    % w = w/sum(w); % normalize weights to sum to 1
    % keyboard
    % Weighted average of real resistivities in frequency range
    % rhoAvg(bb,1) = sum(rho.*w)/sum(w);
    
    % Propagated uncertainty of weighted average
    % rhoAvg(bb,2) = sqrt(sum((w.*drho).^2))/sum(w);
    % rhoAvg(bb,2) = 1/sqrt(sum(w));

    % Kish's design effect variance
    % wbar = sum(w)/N;        
    % w2bar = sum(w.^2)/N;
    % unwVar = (sum(drho.^2))/N;  % unweighted variance
    % wVar = unwVar*w2bar/wbar^2;             % approximate weighted variance
    % rhoAvg(bb,2) = sqrt(wVar);

    % Unbiased standard error
    % nunbias = w2bar*N/((wbar*N)^2-w2bar*N);
    % wVar = sum(w.*drho.^2)/sum(w);
    % rhoAvg(bb,2) = sqrt(wVar*nunbias);

    % % Unweighted average
    % rhoAvg(bb,1) = mean(rho);

    % Propagated error
    % rhoAvg(bb,2) = sqrt(sum(drho.^2))/length(iFind);
    
%     rhoAvg(bb,2) = std(drho)/sqrt(N);
    
end

%% Get rescaled data
if exist('bs_moments_rescaled','var')

    % for cc = 1:length(rescaledOut)
    % 
    %         % Spectral fit temperature
    %         pdT = fitdist(rescaledOut(cc).bootstat(:,1)/t,'normal');
    %         temp(cc,1)=pdT.mu; 
    %         temp(cc,2)=pdT.sigma;
    % 
    %         % Spectral fit gamma
    %         pdG = fitdist(rescaledOut(cc).bootstat(:,2),'normal');
    %         gamma(cc,1)=pdG.mu; 
    %         gamma(cc,2)=pdG.sigma;
    % 
    %         % Spectral fit trap frequency
    %         pdf = fitdist(rescaledOut(cc).bootstat(:,3),'normal');
    %         trap(cc,1)=pdf.mu; 
    %         trap(cc,2)=pdf.sigma;
    % 
    %         % Create parameter arrays
    %         P = [t*temp(cc,1) gamma(cc,1) trap(cc,1)];
    % 
    %         lineshape = conductivity_eval2(FREQ_HZ,P,DEPTH_ER);
    % 
    %         I = find(diff(sign(imag(lineshape))));                                      % index of zero crossings
    %         m = (imag(lineshape(I+1))-imag(lineshape(I)))./(FREQ_HZ(I+1)-FREQ_HZ(I));   % slope
    %         fStarRS(cc,1) = -imag(lineshape(I))./m+FREQ_HZ(I);                                   
    % end  

    %% Average rescaled resistivity
    for bb = 1:length(bs_moments_rescaled)
        
        [M,IMAX] = max(real(bs_moments_rescaled(bb).sigma));
        fStar(bb,1) = bs_moments_rescaled(bb).Frequency_Hz(IMAX);

        % Find frequencies within X% of fStar
        iFind = find(abs(1-bs_moments_rescaled(bb).Frequency_Hz/fStar(bb)) < freqThresh);

        % Real resistivities in frequency range
        rho = real(bs_moments_rescaled(bb).rho(iFind));

        % Uncertainties
        drho = abs(real(bs_moments_rescaled(bb).rhoErr(iFind)/2));

        % Weights
        w = (1./drho).^2;
        % w = w/sum(w); % normalize weights to sum to 1

        % % Weighted average of real resistivities in frequency range
        % rhoAvgRS(bb,1) = sum(rho.*w)/sum(w);
        % 
        % % Propagated uncertainty of weighted average
        % rhoAvgRS(bb,2) = sqrt(sum((w.*drho).^2))/sum(w);

        % % % Unweighted average
        rhoAvgRS(bb,1) = mean(rho);
        % % 
        % % Propagated error
        rhoAvgRS(bb,2) = sqrt(sum(drho.^2))/length(iFind);

    end
end
if exist('bs_moments_Gibbs_rescaled','var')
    %% Average rescaled resistivity
    for cc = 1:length(bs_moments_Gibbs_rescaled)
        
        [M,IMAX] = max(real(bs_moments_Gibbs_rescaled(cc).sigma));
        fStar(cc,1) = bs_moments_Gibbs_rescaled(cc).Frequency_Hz(IMAX);

        % Find frequencies within X% of fStar
        iFind = find(abs(1-bs_moments_Gibbs_rescaled(cc).Frequency_Hz/fStar(cc)) < freqThresh);
        % iFind = find(bs_moments(cc).Frequency_Hz < 100 & bs_moments(cc).Frequency_Hz > 30);
        
        % Real resistivities in frequency range
        rho = real(bs_moments_Gibbs_rescaled(cc).rho(iFind));
        
        % Uncertainties
        drho = abs(real(bs_moments_Gibbs_rescaled(cc).rhoErr(iFind)/2));
        
        % Weights
        w = (1./drho).^2;
        % w = w/sum(w); % normalize weights to sum to 1
    
        % % Weighted average of real resistivities in frequency range
        % rhoAvgGRS(cc,1) = sum(rho.*w)/sum(w);
        % 
        % % Propagated uncertainty of weighted average
        % rhoAvgGRS(cc,2) = sqrt(sum((w.*drho).^2))/sum(w);
    
        % % % Unweighted average
        rhoAvgGRS(cc,1) = mean(rho);
        % % 
        % % Propagated error
        rhoAvgGRS(cc,2) = sqrt(sum(drho.^2))/length(iFind);

    end
end
end

%Define weighted sum function that renormalizes for bootstrap fitting
function [wsum] = rhoAvgfn(wrho_list)
        wsum= sum([wrho_list.rho].*[wrho_list.w])/sum([wrho_list.w]);
        % keyboard
end