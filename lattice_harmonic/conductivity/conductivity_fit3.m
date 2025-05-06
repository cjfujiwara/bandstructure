function output = conductivity_fit3(freq,sigma,sigma_err,LUT,P0)

freq=freq(:);
sigma=sigma(:);
sigma_err=sigma_err(:);

if nargin == 3
    LUT = makeLHO_LUT;
end
N = 101;  

if nargin == 4
    P0 = [700 200 60.65];
end


%% Define Cost Function
 df = LUT.TRAP_Hz(2)-LUT.TRAP_Hz(1);
   
    function yy=error_function(P)
        % P(1) : TEMPERATURE    [Hz]
        % P(2) : GAMMA          [1/s]
        % P(3) : TRAP FREQUENCY [Hz]
        
        % Rename fit parameters 
        
        T       = P(1);
        G       = P(2);
        F       = P(3);

        ind=find(F<LUT.TRAP_Hz,1);
        d1 = F-LUT.TRAP_Hz(ind-1);
        d2 = F-LUT.TRAP_Hz(ind);
        a = 1-abs(d1)/df;
        b = 1-abs(d2)/df;

        ENG = a*LUT.ENG_Hz(:,(ind-1))+b*LUT.ENG_Hz(:,(ind));
        D2 = a*abs(LUT.D2_aL(:,:,ind-1)).^2+b*abs(LUT.D2_aL(:,:,ind)).^2;

        ENG = ENG-ENG(1);
        D2 = D2(1:N,1:N);
        ENG = ENG(1:N);

        % Meshgrid data
        
        [EE1,EE2] = meshgrid(ENG,ENG);
        dEE         = EE1-EE2;       
        % Partition Function
        Z = sum(exp(-ENG/T),'all');
        
        function y = foo(f)
            A = -1i*f*((exp(-EE1/T)-exp(-EE2/T))/Z).*D2./((f-dEE)+1i*G/2/(2*pi));
            y = sum(A,'all');
        end        
         sigma_fit = arrayfun(@(f) foo(f),freq);         
         yy=[(real(sigma_fit)-real(sigma))./real(sigma_err);
         (imag(sigma_fit)-imag(sigma))./imag(sigma_err)]; 
    end


    function [f0,rho0,rhoinf]=getResitivity(P)
      % P(1) : TEMPERATURE    [Hz]
        % P(2) : GAMMA          [1/s]
        % P(3) : TRAP FREQUENCY [Hz]
        
        % Rename fit parameters 
        
        T       = P(1);
        G       = P(2);
        F       = P(3);

        ind=find(F<LUT.TRAP_Hz,1);
        d1 = F-LUT.TRAP_Hz(ind-1);
        d2 = F-LUT.TRAP_Hz(ind);
        a = 1-abs(d1)/df;
        b = 1-abs(d2)/df;

        ENG = a*LUT.ENG_Hz(:,(ind-1))+b*LUT.ENG_Hz(:,(ind));
        D2 = a*abs(LUT.D2_aL(:,:,ind-1)).^2+b*abs(LUT.D2_aL(:,:,ind)).^2;

        ENG = ENG-ENG(1);
        D2 = D2(1:N,1:N);
        ENG = ENG(1:N);

        % Meshgrid data
        
        [EE1,EE2] = meshgrid(ENG,ENG);
        dEE         = EE1-EE2;       
        % Partition Function
        Z = sum(exp(-ENG/T),'all');
        
        function y = foo(f)
            A = -1i*f*((exp(-EE1/T)-exp(-EE2/T))/Z).*D2./((f-dEE)+1i*G/2/(2*pi));
            y = sum(A,'all');
        end     

        % Find Resitivity at zero imag cond
        fvec = linspace(F/4,2*F,100);
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
% options.StepTolerance = 0.1;
options.Display='off';
options.FiniteDifferenceStepSize=0.01;
% options = optimoptions('lsqnonlin', 'FiniteDifferenceStepSize', mesh_step_size);

[fout, resnorm, residual, exitflag, output0,...
    lambda, jacobian] = lsqnonlin(@(P) error_function(P), P0,[],[],options);
conf = nlparci(fout,residual,'jacobian',jacobian);

[f0,rho0, rhoinf]=getResitivity(fout);

%% Create Output

output = struct;
% output.freq_fit = vals.freq_fit;
% output.sigma_fit = vals.sigma_fit;
output.fout = fout;
output.conf=conf;
output.f0 = f0;
output.rho0=rho0;
output.rhoinf=rhoinf;

end

