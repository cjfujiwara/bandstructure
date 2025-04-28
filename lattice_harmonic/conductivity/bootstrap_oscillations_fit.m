function output=bootstrap_oscillations_fit(x,y1,y2,f)

% Amplitude Guess
guess_Amplitude = 0.5*(max(y1)-min(y1));

% Center Guess
guess_x0 = median(y1);

% Velocity Guess
guess_v0 = 0;

% Phase Guess
phiVec = linspace(0,-2*pi,100);
phi_sse = zeros(length(phiVec),1);
for rr=1:length(phiVec)
    phi_sse(rr) = sum((guess_Amplitude*sin(2*pi*f*x + phiVec(rr)) + guess_x0 - y1).^2);
end
[~,ii] = min(phi_sse);
guess_phi = phiVec(ii);

S_guess = -guess_Amplitude*cos(guess_phi);
C_guess = -guess_Amplitude*sin(guess_phi);
Tc = (max(x)+min(x))*0.5;

oscillations_wrapper = @(P,t) ...
    -P(1)*sin(2*pi*f*t) ...
    -P(2)*cos(2*pi*f*t) + ...
    P(3) + P(4)*(t-Tc);

P_guess = [S_guess C_guess guess_x0 guess_v0];
yG = oscillations_wrapper(P_guess,x);
    
%% Normal Fitting
        options = optimset('Display','off');    
        x = x(:);
y1 = y1(:);
data = [x y1];

[fout,resnorm,residual,exitflag,output0,lambda,jacobian] = ...
    lsqcurvefit(oscillations_wrapper, P_guess, data(:,1), data(:,2), [], [], options);
conf = nlparci(fout,residual,'jacobian',jacobian);
% SS_res = resnorm;
% SS_tot = sum((real(z)-mean(real(z))).^2+(imag(z)-mean(imag(z))).^2);
% R2 = 1 - SS_res/SS_tot;

S_val = fout(1);
S_err = (conf(1,2)-conf(1,1))/2;

C_val = fout(2);
C_err = (conf(2,2)-conf(2,1))/2;

x0_val = fout(3);
x0_err = (conf(3,2)-conf(3,1))/2;

v0_val = fout(4);
v0_err = (conf(4,2)-conf(4,1))/2;

P_fit=[S_val C_val x0_val v0_val];
P_err=[S_err C_err x0_err v0_err];

%% Bootstrap Fitting
    function fittedParams=fitModel(data)
        options = optimset('Display','off');    
        fittedParams = lsqcurvefit(oscillations_wrapper, P_guess, data(:,1), data(:,2), [], [], options);
    end

nBootstraps = 1e3;
% Apply bootstrap
[bootstat, bootsam] = bootstrp(nBootstraps, @fitModel, data);

figure;
clf
subplot(241)
histfit(bootstat(:,1));
pdS = fitdist(bootstat(:,1),'Normal');
xlabel('S (um)');
ylabel('occurences')
str=['boot normal : ' num2str(round(pdS.mu,2)) '\pm' num2str(round(2*pdS.sigma,2)) ...
    newline ...
    'fit 95 conf : ' num2str(round(S_val,2)) '\pm' num2str(round(S_err,2))];
text(.01,.99,str,'units','normalized','verticalalignment','top','backgroundcolor',[1 1 1 .5]);

subplot(242)
histfit(bootstat(:,2));
pdC = fitdist(bootstat(:,2),'Normal');
xlabel('C (um)');
ylabel('occurences')
str=['boot normal : ' num2str(round(pdC.mu,2)) '\pm' num2str(round(2*pdC.sigma,2)) ...
    newline ...
    'fit 95 conf : ' num2str(round(C_val,2)) '\pm' num2str(round(C_err,2))];
text(.01,.99,str,'units','normalized','verticalalignment','top','backgroundcolor',[1 1 1 .5]);

subplot(243)
histfit(bootstat(:,3));
xlabel('x0 (um)');
ylabel('occurences')

subplot(244)
histfit(bootstat(:,4));
xlabel('v0 (um/ms)');
ylabel('occurences')

subplot(2,4,[5 6 7])
co=get(gca,'colororder');
tt=linspace(min(x),max(x),100);
plot(tt,oscillations_wrapper(P_fit,tt),'r-');
hold on
plot(x,y1,'o','markerfacecolor',co(1,:),'color','k');
hold on
xlabel('total time (ms)');
ylabel('position (um)')

subplot(2,4,8)
plot(y2,y1,'o');
xlabel('y1')
ylabel('y2');
axis equal tight
%% Create Ouputs

output = struct;
output.x = x;
output.y = y1;
output.FitParam = P_fit;
output.FitErr = P_err;
output.BootStat = bootstat;
output.BootSam = bootsam;
end

