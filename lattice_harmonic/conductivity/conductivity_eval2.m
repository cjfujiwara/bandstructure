function z = conductivity_eval2(FREQ_HZ,P,DEPTH_ER)
% Author : CJ Fujiwara
%
% This function evaluates the conductivity assuming lineare response and a
% broadening parameter Gamma.

if nargin == 2
    DEPTH_ER=2.5;
end
%% Convert Fit Param into each value
TEMP_HZ = P(1);
GAMMA_INV_SEC = P(2);
TRAP_HZ     = P(3);

%% Settings

Nsites = 601;
TunnelOrder=15;
HarmonicBands=1;

%% Calculate Lattice Properties

if isstruct(DEPTH_ER)
    lattice = DEPTH_ER;
else

    lattice                     = constants;
    lattice.depth               = DEPTH_ER; 
    lattice.WannierBands        = [1];
    
    lattice.numStates           = 101;       % must be odd
    lattice.numK                = 301;      % must be odd    
    
    wannier_opts                = struct;
    wannier_opts.Bands          = [1];
    
    
    lattice = calculateBandStructure(lattice);      % calculate band structure
    lattice = calculateTunneling(lattice);          % calculate tunneling elements
    lattice = wannier(lattice,wannier_opts);        % Calculate wannier function
    lattice = calculateWannierMoments(lattice);     % Dipole matrix elements in wannier basis
end
%%

 function sigma_fit=sigma_eval(P)
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
        sigma_fit = arrayfun(@(f) foo(f),FREQ_HZ);         
 end


z = sigma_eval(P);

end

