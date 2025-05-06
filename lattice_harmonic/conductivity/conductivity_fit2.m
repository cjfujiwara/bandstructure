function output = conductivity_fit2(lattice,freq,sigma,sigma_err,input)

% freq      : frequency data
% sigma     : complex conductivity data

if nargin ==2
    sigma_err = zeros(length(freq),1);
end

%% Settings

% Number of eigenstates to include in fit
N = 101;  
% Numerical Settings
Nsites          = 301;
TunnelOrder     = 11;
HarmonicBands   = 1;    

%% Process Data

% Make sure they are a colummn vector
freq = freq(:);     
sigma = sigma(:);
sigma_err = sigma_err(:);

%% Construct Initial Guess

P0 = [700 200 57];
if nargin==5
   P0=input.fout; 
end

%% Define Cost Function
        
    function yy=error_function(P)
        % P(1) : TEMPERATURE   [Hz]
        % P(2) : GAMMA         [1/s]
        % P(3) : TRAP FREQUENC [Hz]
        
        % Rename fit parameters 
        T       = P(1);
        G       = P(2);
        omega   = 2*pi*P(3); 

        % disp(num2str(P(3)))
        
        
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
        
         sigma_fit = arrayfun(@(f) foo(f),freq);
         
         yy=[(real(sigma_fit)-real(sigma))./real(sigma_err);
         (imag(sigma_fit)-imag(sigma))./imag(sigma_err)];        
    end

    function [rho0,rhoinf]=getVals(P)
        % Rename fit parameters 
        T       = P(1);
        G       = P(2);
        omega   = 2*pi*P(3);        
        F=P(3);
        
        
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
        fvec = linspace(F/4,2*F,50);
        sigma_fit = arrayfun(@(f) foo(f),fvec);         
        ig=find(sign(diff(sign(imag(sigma_fit))))==1,1);         
        f0=fvec(ig);         
        f0=fzero(@(f) imag(foo(f)),f0);  
        rho0=real(1/foo(f0));      

        % find high freq rho
        rhoinf = real(1/foo(1e3));     
    end

%% Fit it
%        options = optimset('Display','on','TolFun',1e-3);    
options = optimoptions('lsqnonlin');
options.FunctionTolerance=1e-9;
options.OptimalityTolerance=1e-9;
options.Display='off';

[fout, resnorm, residual, exitflag, output0,...
    lambda, jacobian] = lsqnonlin(@(P) error_function(P), P0,[],[],options);
conf = nlparci(fout,residual,'jacobian',jacobian);

[rho0,rhoinf]=getVals(fout);

%% Create Output

output = struct;
output.fout = fout;
output.conf=conf;
output.rho0=rho0;
output.rhoinf=rhoinf;

end

