function output = conductivity_fit2(lattice,FREQ,SIGMA,SIGMA_ERR,input)

% freq      : frequency data
% sigma     : complex conductivity data
t1=now;
if nargin ==2
    SIGMA_ERR = zeros(length(FREQ),1);
end

%% Settings

% Numerical Settings
Nsites          = 301;      % number of lattice sites to consider
TunnelOrder     = 11;       % tunneling order to consider
HarmonicBands   = 1;        % which bands to use (always 1 for now)
N               = 101;      % Number of eigenstates to include in fit

%% Process Data

% Make sure they are a column vector
FREQ        = FREQ(:);     
SIGMA       = SIGMA(:);
SIGMA_ERR   = SIGMA_ERR(:);

%% Construct Initial Guess
% WE NEED TO FIT A DRUDE/LORENTZIAN AND THEN MAKE A SMART GUESS
% SUM RULE GETS US TEMP
% FWHM GETS GAMMA
if nargin==5
   P0=input.fout; 
else
    fprintf('calculating initial guess')
    drude = conductivity_fit_drude(FREQ,SIGMA,SIGMA_ERR);    
    % Band mass at k=0;
    m0          = lattice.BandMassGamma(1); % in units of bare mass
    amu         = 1.66054e-27;      % [kg] atomic mass unit
    m           = 40*amu;           % [kg] potassium-40 mass
    aL          = 532e-9;           % [m] lattice spacing
    hbar        = 1.05457182e-34;   % [Js] reduce planck's constant
    S_DRUDE     = drude.SumRule;    
    t           = lattice.Tunneling(1,1)*lattice.fr;

    % Match the peak conductivity at this Gamma to make initial guess
    tVec        = [0.1:.1:2 2.5:.5:6];
    gamma_drude = abs(drude.fout(2));
    f0_drude    = drude.fout(3)*sqrt(m0);


    f0_drude    = min([f0_drude 61]);
    s0_DRUDE    = drude.SIGMA_PEAK;
    FREQ_PEAK   = drude.fout(3);
    s0_TDPT = zeros(length(tVec),1);
    for jj = 1:length(tVec)
        P0 = [tVec(jj)*t gamma_drude f0_drude];
        s0_TDPT(jj)  = real(conductivity_eval2(FREQ_PEAK, P0,lattice));
    end
    Tg = interp1(s0_TDPT,tVec,s0_DRUDE);
    P0 = [Tg*t gamma_drude f0_drude];    
    fprintf('done')

    % uncomment for debugging
    % FREQ_THEORY = [0:2:200 250:50:500 2e3];
    % plot(drude.FREQ_THEORY,real(drude.SIGMA_THEORY));hold on;
    % sigma_TDPT  = conductivity_eval2(FREQ_THEORY, P0,lattice);
    % plot(FREQ_THEORY,real(sigma_TDPT));
end




%% Make Frequency matrix
% N x N x n (where n = length of frequencies)

FREQ_MAT = zeros(N,N,length(FREQ));
for jj=1:length(FREQ)
    a = FREQ(jj);
    FREQ_MAT(:,:,jj) = a(ones(N,N));
end

%% Define Cost Function

    function yy=error_function(P)
        ta=now;
        % P(1) : TEMPERATURE   [Hz]
        % P(2) : GAMMA         [1/s]
        % P(3) : TRAP FREQUENCY [Hz]
        
        % Rename fit parameters 
        T       = P(1);
        G       = P(2);
        omega   = 2*pi*P(3); 
        
        % Solve Eigenvalue Problem
        opts=struct;               
        opts.NumSites = Nsites;
        opts.MaxTunnelingOrder = TunnelOrder;
        opts.HarmonicBands = HarmonicBands;
        opts.omega = omega;
        opts.Omega = 0.5*lattice.m*omega^2*(lattice.lambda/2)^2/lattice.h;
        [~,lho] = calculateLHOSpectrum_sband(lattice,opts);
        
        % Load Eigenvalues and Dipole Operator
        N           = 101;
        d2          = abs(lho.DipoleOperator).^2;
        eng         = lho.EigenValues;
        eng         = eng-eng(1);    
        eng         = eng(1:N);
        d2          = d2(1:N,1:N);

        % Meshgrid data
        [myEE1,myEE2] = meshgrid(eng,eng);
        mydEE         = myEE1-myEE2;   
    
        % Partition Function
        Z = sum(exp(-eng/T),'all');

        % OLD WAY
        % sigma[freq_drive] for a single frequency
        % function y = foo(f)
        %     A = -1i*f*((exp(-myEE1/T)-exp(-myEE2/T))/Z).*d2./((f-mydEE)+1i*G/2/(2*pi));
        %     y = sum(A,'all');
        % end        
        % % sigma_fit = arrayfun(@(f) foo(f),FREQ); % evaluate sigma for all drive freqs
        % 
        % sigma_fit = zeros(length(FREQ),1);
        % for cc=1:length(FREQ)
        %     sigma_fit(cc) = foo(FREQ(cc));
        % end
        
        % Matrix way
        mydEE3 = repmat(mydEE,[1 1 length(FREQ)]);
        myEE13  = repmat(myEE1,[1 1 length(FREQ)]);
        myEE23  = repmat(myEE2,[1 1 length(FREQ)]);
        D23     = repmat(d2,[1 1 length(FREQ)]);
        A3 = -1i*FREQ_MAT.*((exp(-myEE13/T)-exp(-myEE23/T))/Z).*D23./((FREQ_MAT-mydEE3)+1i*G/2/(2*pi));
        sigma_fit = sum(A3,[1 2]);
        sigma_fit = sigma_fit(:);

         yy=[(real(sigma_fit)-real(SIGMA))./real(SIGMA_ERR);
         (imag(sigma_fit)-imag(SIGMA))./imag(SIGMA_ERR)]; 
        tb=now;
        % disp((tb-ta)*24*60*60)

        
    end

    function [rho0,rhoinf]=getVals(P)
        % Rename fit parameters 
        T       = P(1);
        G       = P(2);
        omega   = 2*pi*P(3);        
        F=P(3);
        ta=now;

        
        % Solve Eigenvalue Problem
        opts=struct;               
        opts.NumSites = Nsites;
        opts.MaxTunnelingOrder = TunnelOrder;
        opts.HarmonicBands = HarmonicBands;
        opts.omega = omega;
        opts.Omega = 0.5*lattice.m*omega^2*(lattice.lambda/2)^2/lattice.h;
        [~,lho] = calculateLHOSpectrum_sband(lattice,opts);
        
        % Load Eigenvalues and Dipole Operator
        N           = 101;
        d2          = abs(lho.DipoleOperator).^2;
        eng         = lho.EigenValues;
        eng         = eng-eng(1);    
        eng         = eng(1:N);
        d2          = d2(1:N,1:N);

        % Meshgrid data
        [myEE1,myEE2] = meshgrid(eng,eng);
        mydEE         = myEE1-myEE2;   
    
        % Partition Function
        Z = sum(exp(-eng/T),'all');
        
        function y = foo(f)
            A = -1i*f*((exp(-myEE1/T)-exp(-myEE2/T))/Z).*d2./((f-mydEE)+1i*G/2/(2*pi));
            y = sum(A,'all');
        end
        
        % Find Resitivity at zero imag cond
        fvec = linspace(F/4,2*F,20);
        sigma_fit = arrayfun(@(f) foo(f),fvec);         
        ig=find(sign(diff(sign(imag(sigma_fit))))==1,1);         
        f0=fvec(ig);     
        if isempty(f0)
            f0 = 60;
        end
        f0=fzero(@(f) imag(foo(f)),f0);
        rho0=real(1/foo(f0));      
        % rho0 = 0.3;

        % find high freq rho
        rhoinf = real(1/foo(1e3)); 
        % rhoinf = 0.3;
        tb=now;
        % disp((tb-ta)*24*60*60)
    end

%% Fit it
%        options = optimset('Display','on','TolFun',1e-3);    
options = optimoptions('lsqnonlin');
options.FunctionTolerance=1e-9;
options.OptimalityTolerance=1e-9;
options.Display='off';
lb = [1000 1 64];[0 0 0];
ub = [3500 inf 64];

[fout, resnorm, residual, exitflag, output0,...
    lambda, jacobian] = lsqnonlin(@(P) error_function(P), P0,lb,ub,options);
conf = nlparci(fout,residual,'jacobian',jacobian);


[rho0,rhoinf]=getVals(fout);
t2=now;

%% Create Output

% disp((t2-t1)*24*60*60)

% keyboard
output = struct;
output.fout = fout;
output.conf=conf;
output.rho0=rho0;
output.rhoinf=rhoinf;

end