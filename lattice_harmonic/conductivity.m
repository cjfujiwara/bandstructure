function fout=conductivity(x,z,src)

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
    omega = 2*pi*[60:1:75];
    
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
    function y = sigma(T,G,f)
        y = arrayfun(@(f) sigma_helper(T,G,f),f);
    end

    function y = sigma_helper(T,G,f)
        Z = sum(exp(-eng/T),'all');
        A = -1i*f*((exp(-EE1/T)-exp(-EE2/T))/Z).*d2./((f-dEE)+1i*G/2/(2*pi));
        y = sum(A,'all');
    end

% Real conductivity
    function y = sigma_real(T,G,f)        
        y = real(sigma(T,G,f));
    end

% Imaginary conductivity
    function y = sigma_imag(T,G,f)
        y = imag(sigma(T,G,f));
    end

%%

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
P = [1100 51];
% P = []
%% Constrained Fit
y = [real(z); imag(z)];
foo = @(P,f) [sigma_real(P(1),P(2),f); sigma_imag(P(1),P(2),f)];

[fout,resnorm,residual,exitflag,output0,lambda,jacobian]=lsqcurvefit(foo,P,x,y);
conf = nlparci(fout,residual,'jacobian',jacobian);

SS_res = resnorm;
SS_tot = sum((real(z)-mean(real(z))).^2+(imag(z)-mean(imag(z))).^2);
R2 = 1 - SS_res/SS_tot;

T = fout(1);
Terr = (conf(1,2)-conf(1,1))/2;

G = fout(2);
Gerr = (conf(2,2)-conf(2,1))/2;

output = struct;
output.T = T;
output.Terr = Terr;
output.G = G;
output.Gerr = Gerr;
output.Rsquared = R2;

disp(output)
%
ft = linspace(0,100,500);
yt = sigma(T,G,ft);

figure(21);
% clf
subplot(121);
plot(ft,real(yt),'r-');
hold on
plot(x,real(z),'ko','markerfacecolor','k');
subplot(122);
plot(ft,imag(yt),'r-');
hold on
plot(x,imag(z),'ko','markerfacecolor','k');
% disp(resnorm);

%% Real Fit

%% Imaginary Fit


end

