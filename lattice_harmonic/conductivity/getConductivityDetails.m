function [outputArg1,outputArg2] = getConductivityDetails(P,LUT)

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

    fvec=linspace(1,300,1e3);        
    sigma_fit = arrayfun(@(f) foo(f),fvec);         
    ig=find(sign(diff(sign(imag(sigma_fit))))==1,1);         
    f0=fvec(ig);         
    f0=fzero(@(f) imag(foo(f)),f0);  

    rho0=real(1/foo(f0));         
    rhoinf = real(1/foo(1e3));        

    vals=struct;
    vals.rho0=rho0;
    vals.rhoinf=rhoinf;
    vals.freq_fit = fvec;
    vals.sigma_fit = sigma_fit;        
    

end

