function fout = conductivity(x,z,src)

y = [real(z); imag(z)];

figNum1 = 1990;
figNum2 = 1991;

%% Drude Fit
drude_wrapper = @(P,f) [real(drude(P(1),P(2),P(3),f)); imag(drude(P(1),P(2),P(3),f))];

function y = drude(amp,f0,G,f)
    w   = 2*pi*f;
    w0  = 2*pi*f0;
    y   = 1i*amp*G*w./((w.^2-w0.^2)+1i*w*G);
end
[amp,ind] = max(real(z));
f0 = x(ind);
P_drude_complex = [amp f0 50];

[fout,resnorm,residual,exitflag,output0,lambda,jacobian]= ...
    lsqcurvefit(drude_wrapper,P_drude_complex,x,y);
conf = nlparci(fout,residual,'jacobian',jacobian);

SS_res = resnorm;
SS_tot = sum((real(z)-mean(real(z))).^2+(imag(z)-mean(imag(z))).^2);
R2 = 1 - SS_res/SS_tot;

amp = fout(1);
amperr = (conf(1,2)-conf(1,1))/2;

f0 = fout(2);
f0err = (conf(2,2)-conf(2,1))/2;

G = fout(3);
Gerr = (conf(3,2)-conf(3,1))/2;

ft = linspace(0,1000,1000);

drude_complex               = struct;
drude_complex.s0            = amp;
drude_complex.s0err         = amperr;
drude_complex.f0            = f0;
drude_complex.f0err         = f0err;
drude_complex.G             = G;
drude_complex.Gerr          = Gerr;
drude_complex.Rsquared      = R2;
drude_complex.FreqFit       = ft;
drude_complex.SigmaFit      = drude(amp,f0,G,ft);

drude_sum = trapz(ft,real(drude(amp,f0,G,ft)));

% keyboard

doCalculateLHO = 1;
if doCalculateLHO
    %% Lattice Properties    
    lattice                 = constants;
    lattice.depth           = [2.5]; 
    lattice.WannierBands    = [1];
    
    lattice.numStates       = 101;       % must be odd
    lattice.numK            = 301;      % must be odd    
    
    wannier_opts            = struct;
    wannier_opts.Bands      = [1];
    
    lattice = calculateBandStructure(lattice);      % calculate band structure
    lattice = calculateTunneling(lattice);          % calculate tunneling elements
    lattice = wannier(lattice,wannier_opts);        % Calculate wannier function
    lattice = calculateWannierMoments(lattice);     % Dipole matrix elements in wannier basis
    
    %% Harmonic Properties
    
    % Numerical Settings
    Nsites = 601;
    TunnelOrder = 21;
    HarmonicBands = 1;    
    % omega = 2*pi*[65:1:75];
        omega = 2*pi*[(f0+2):1:(f0+12)];

    %% Calculate Lattice + HO Properties
    
    HO_opts = struct;
    HO_opts.NumSites = Nsites;
    HO_opts.MaxTunnelingOrder = TunnelOrder;
    HO_opts.HarmonicBands = HarmonicBands;
    for nn=1:length(omega)
        HO_opts.omega = omega(nn);
        HO_opts.Omega = 0.5*lattice.m*omega(nn)^2*(lattice.lambda/2)^2/lattice.h;
        [lattice,LHO(nn)] = calculateLHOSpectrum_sband(lattice,HO_opts);
    end    
    src=LHO(5);
end

%% Conductivity Functions

% Total conductivity
    function y = sigma_func(T,G,f)
        y = arrayfun(@(f) sigma_helper(T,G,f),f);
    end

    function y = sigma_helper(T,G,f)
        Z = sum(exp(-eng/T),'all');
        A = -1i*f*((exp(-EE1/T)-exp(-EE2/T))/Z).*d2./((f-dEE)+1i*G/2/(2*pi));
        y = sum(A,'all');
    end

% Real conductivity
    function y = sigma_real(T,G,f)        
        y = real(sigma_func(T,G,f));
    end

% Imaginary conductivity
    function y = sigma_imag(T,G,f)
        y = imag(sigma_func(T,G,f));
    end

%%
output = struct;
hF1 = figure(figNum1);
hF1.Color='w';
clf


ax1 = subplot(2,2,1,'parent',hF1);
ax2 = subplot(2,2,2,'parent',hF1);
ax3 = subplot(2,2,3,'parent',hF1);
ax4 = subplot(2,2,4,'parent',hF1);
cc=jet(length(LHO));

for jj=1:length(LHO)
    src = LHO(jj);
    
    % Number of eigenstates to include in fit
    N = 101;
    
    % Make into column vector
    x=x(:);
    z=z(:);
    
    % Load Eigenvalues and Dipole Operator
    d2 = abs(src.DipoleOperator).^2;
    eng = src.EigenValues;
    eng = eng-eng(1);
    
    % Reduce vector space
    eng=eng(1:N);
    d2=d2(1:N,1:N);
    
    % Meshgrid data
    [EE1,EE2] = meshgrid(eng,eng);
    dEE = EE1-EE2;
    
    P = [850 200];
    % P = [1100 51];
    
    %% Constrained Fit
    y = [real(z); imag(z)];
    drude_wrapper = @(P,f) [sigma_real(P(1),P(2),f); sigma_imag(P(1),P(2),f)];
    
    [fout,resnorm,residual,exitflag,output0,lambda,jacobian]=lsqcurvefit(drude_wrapper,P,x,y);
    conf = nlparci(fout,residual,'jacobian',jacobian);
    
    SS_res = resnorm;
    SS_tot = sum((real(z)-mean(real(z))).^2+(imag(z)-mean(imag(z))).^2);
    R2 = 1 - SS_res/SS_tot;
    
    T = fout(1);
    Terr = (conf(1,2)-conf(1,1))/2;
    
    G = fout(2);
    Gerr = (conf(2,2)-conf(2,1))/2;
    output(jj).trap_freq = src.omega/(2*pi);
    output(jj).T = T;
    output(jj).Terr = Terr;
    output(jj).G = G;
    output(jj).Gerr = Gerr;
    output(jj).Rsquared = R2;
    
    disp(output)
    %
    % ft = linspace(0,100,500);
    yt = sigma_func(T,G,ft);  
    plot(ft,real(yt),'-','color',cc(jj,:),'parent',ax1);
    hold(ax1,'on')
    plot(ft,imag(yt),'-','color',cc(jj,:),'parent',ax2);
    hold(ax2,'on')
    % disp(resnorm);
end

plot(x,real(z),'ko','markerfacecolor','k','parent',ax1);
plot(x,imag(z),'ko','markerfacecolor','k','parent',ax2);

subplot(2,2,3,'parent',hF1);
plot(omega/(2*pi),[output.Rsquared],'o','parent',ax3)
xlabel(ax3,'trap freq (Hz)')
ylabel(ax3,'R squared')

subplot(2,2,4,'parent',hF1);
yyaxis left
errorbar(omega/(2*pi),[output.G],[output.Gerr],'o','parent',ax4)
xlabel(ax4,'trap freq (Hz)')
ylabel(ax4,'\Gamma (1/s)')
yyaxis right
errorbar(omega/(2*pi),[output.T]/563,[output.Terr]/563,'o','parent',ax4)
ylabel('temp (t)')
ylim([0 3])

%% Summary Figure
[val,ind]=max([output.Rsquared]);
output_best = output(ind);
T=output_best.T;
Terr=output_best.Terr;
G=output_best.G;
Gerr=output_best.Gerr;
f0_best = output_best.trap_freq;

t=LHO(1).Tunneling(1);
%%
hF2 = figure(figNum2);
hF2.Color='w';
clf(hF2);
axes('parent',hF2)
co=get(gca,'colororder');
% axa = subplot(2,1,1,'parent',hF2);

pDrude_R=plot(ft,real(drude(amp,f0,G,ft)),'--','color',co(1,:));
hold on
pDrude_I=plot(ft,imag(drude(amp,f0,G,ft)),'--','color',co(2,:));
ylabel('$\sigma/\sigma_0$','interpreter','latex')
xlabel('drive frequency (Hz)')
set(gca,'box','on','fontname','times')

pTDPT_R=plot(ft,real(sigma_func(T,G,ft)),'-','color',co(1,:),'linewidth',1);
hold on
pTDPT_I=plot(ft,imag(sigma_func(T,G,ft)),'-','color',co(2,:),'linewidth',1);

pData_R=plot(x,real(z),'o','markerfacecolor',co(1,:),'color',co(1,:)*.5);
pData_I=plot(x,imag(z),'s','markerfacecolor',co(2,:),'color',co(2,:)*.5);

str1=['Re[drude]'];
str2=['Im[drude]'];

str3=['Re[tdpt]'];
str4=['Im[tdpt]'];

legStr={str1, str2, str3, str4};
legend([pDrude_R pDrude_I pTDPT_R pTDPT_I],legStr)

str_drude = ['drude $(' ...
    'f_0:' num2str(round(drude_complex.f0,1)) '~\mathrm{Hz}' ...
    ',A:' num2str(round(drude_complex.s0,1))  ...    
    ',\Gamma:' num2str(round(drude_complex.G,1)) '\mathrm{s}^{-1})$'];

str_tdpt = ['tdpt $(f_0:' num2str(round(f0_best,1)) '~\mathrm{Hz}' ...
    ',T:' num2str(round(T/t,2)) 't'...
    ',\Gamma:' num2str(round(G,1)) '\mathrm{s}^{-1})$'];

str_tdpt_2 = ['tdpt err $\pm' num2str(round(Terr/t,1)) 't'...
    ',\pm' num2str(round(Gerr,0)) '\mathrm{s}^{-1})$'];

str = [str_drude newline str_tdpt newline str_tdpt_2];

text(.99,.01,str,'units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom','interpreter','latex')
xlim([0 150]);
% keyboard

s='Mar16spectrum198p5';
%% Real Fit
%%
hF3 = figure(figNum2+1);
hF3.Color='w';
clf(hF3);
axes('parent',hF3)
co=get(gca,'colororder');
% axa = subplot(2,1,1,'parent',hF2);

pDrude_R=plot(ft,real(1./drude(amp,f0,G,ft)),'--','color',co(1,:));
hold on
pDrude_I=plot(ft,imag(1./drude(amp,f0,G,ft)),'--','color',co(2,:));
ylabel('$\rho$','interpreter','latex')
xlabel('drive frequency (Hz)')
set(gca,'box','on','fontname','times')

pTDPT_R=plot(ft,real(1./sigma_func(T,G,ft)),'-','color',co(1,:),'linewidth',1);
hold on
pTDPT_I=plot(ft,imag(1./sigma_func(T,G,ft)),'-','color',co(2,:),'linewidth',1);

pData_R=plot(x,real(1./z),'o','markerfacecolor',co(1,:),'color',co(1,:)*.5);
pData_I=plot(x,imag(1./z),'s','markerfacecolor',co(2,:),'color',co(2,:)*.5);

str1=['Re[drude]'];
str2=['Im[drude]'];

str3=['Re[tdpt]'];
str4=['Im[tdpt]'];

legStr={str1, str2, str3, str4};
legend([pDrude_R pDrude_I pTDPT_R pTDPT_I],legStr)

str_drude = ['drude $(' ...
    'f_0:' num2str(round(drude_complex.f0,1)) '~\mathrm{Hz}' ...
    ',A:' num2str(round(drude_complex.s0,1))  ...    
    ',\Gamma:' num2str(round(drude_complex.G,1)) '\mathrm{s}^{-1})$'];

str_tdpt = ['tdpt $(f_0:' num2str(round(f0_best,1)) '~\mathrm{Hz}' ...
    ',T:' num2str(round(T/t,2)) 't'...
    ',\Gamma:' num2str(round(G,1)) '\mathrm{s}^{-1})$'];

str_tdpt_2 = ['tdpt err $\pm' num2str(round(Terr/t,1)) 't'...
    ',\pm' num2str(round(Gerr,0)) '\mathrm{s}^{-1})$'];

str = [str_drude newline str_tdpt newline str_tdpt_2];

text(.99,.01,str,'units','normalized','horizontalalignment','right',...
    'verticalalignment','bottom','interpreter','latex')
xlim([20 100]);
% keyboard

s='Mar16spectrum198p5';
%% Imaginary Fit


figure(7);
clf
set(gcf,'Color','w');

i1 = trapz(ft,real(sigma_func(563,250,ft)));
i2 = trapz(ft,real(sigma_func(563*1.5,250,ft)));
i3 = trapz(ft,real(sigma_func(563*2.0,250,ft)));
i4 = trapz(ft,real(sigma_func(563*2.5,250,ft)));

subplot(121);
plot(ft,real(1./sigma_func(563,250,ft))); hold on
plot(ft,real(1./sigma_func(563*1.5,250,ft))); hold on
plot(ft,real(1./sigma_func(563*2.0,250,ft))); hold on
plot(ft,real(1./sigma_func(563*2.5,250,ft))); hold on
ylabel('Re(\sigma_0/\sigma)')
legend({'T=1.0t','T=1.5t','T=2.0t','T=2.5t'})
xlim([0 100])

subplot(122);
plot(ft,real(1./sigma_func(563*3,50,ft))); hold on
plot(ft,real(1./sigma_func(563*3,100,ft))); hold on
plot(ft,real(1./sigma_func(563*3,150,ft))); hold on
plot(ft,real(1./sigma_func(563*3,200,ft))); hold on
ylabel('Re(\sigma_0/\sigma)')
legend({'\Gamma=50/s','T=100/s','T=150/s','T=200/s'})

%%

keyboard

end

