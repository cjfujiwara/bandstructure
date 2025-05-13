function output = conductivity_fit_drude(FREQ,SIGMA,SIGMA_ERR)
%% Formatting


% Make sure they are a colummn vector
FREQ        = FREQ(:);     
SIGMA       = SIGMA(:);
SIGMA_ERR   = SIGMA_ERR(:);


%% Conductivity Function
    function y = sigma_real_fit(A,G,omega0,omega)
        y = A*((G.*omega).^2./((omega.^2-omega0.^2).^2+(omega.*G).^2));
    end
    function y = sigma_imag_fit(A,G,omega0,omega)
        y =  A*G.*((omega.^2-omega0.^2).*omega./((omega.^2-omega0.^2).^2+(omega*G).^2));     
    end

%% Make Guess

[Ag,ind]=max(real(SIGMA));
f0g= FREQ(ind);

Gg = 0.5*2*pi*2.3548*sqrt(sum((FREQ-f0g).^2.*real(SIGMA)/sum(real(SIGMA))));

P0 = [Ag Gg f0g];

%% Fit it
options = optimoptions('lsqnonlin');
options.FunctionTolerance=1e-9;
options.OptimalityTolerance=1e-9;
options.Display='off';

[fout, resnorm, residual, exitflag, output0,...
    lambda, jacobian] = lsqnonlin(@(P) error_function(P), P0,[],[],options);
conf = nlparci(fout,residual,'jacobian',jacobian);

%% Cost Function
    function yy=error_function(P)
        % P(1) : AMPLITUDE      [sigma_0]
        % P(2) : GAMMA          [1/s]
        % P(3) : TRAP FREQUENCY [Hz]
        
        % Rename fit parameters 
        A       = P(1);
        G       = P(2);
        omega0  = 2*pi*P(3);
        omega   = 2*pi*FREQ;   
         
         yy = [(sigma_real_fit(A,G,omega0,omega)-real(SIGMA))./real(SIGMA_ERR);
            (sigma_imag_fit(A,G,omega0,omega)-imag(SIGMA))./imag(SIGMA_ERR)];       
    end

%% Calculate Sum Rule
ff = linspace(0,2e3,1e4);

Sr = (2/pi)*trapz(2*pi*ff,sigma_real_fit(fout(1),fout(2),2*pi*fout(3),2*pi*ff));

s = sigma_real_fit(fout(1),fout(2),2*pi*fout(3),2*pi*ff)+...
    1i*sigma_imag_fit(fout(1),fout(2),2*pi*fout(3),2*pi*ff);

%% Create Output

output = struct;
output.fout = fout;
output.conf=conf;
output.SumRule = Sr;
output.FREQ_THEORY = ff;
output.SIGMA_THEORY = s;
output.SIGMA_PEAK = max(real(s));
end

