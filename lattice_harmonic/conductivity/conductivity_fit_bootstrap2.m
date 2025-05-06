function out=conductivity_fit_bootstrap2(freq,sigma,sigma_err,LUT)


data=[freq sigma sigma_err];

%% Calculate Lattice Properties


%% Fit it once

normal_fit = conductivity_fit3(freq,sigma,sigma_err,LUT);


%%
n=0;
     function fittedParams=fitModel(data)
         freq = data(:,1);
        sigma = data(:,2);
        sigma_err = data(:,3);         
        output = conductivity_fit3(freq,sigma,sigma_err,LUT);
        fittedParams= [output.fout output.rho0 output.rhoinf]; 
        n=n+1;
        % disp([num2str(n) ' ' num2str(t2) ' sec.']);
     end


nBootstraps = 1000;
% Apply bootstrap
[bootstat, bootsam] = bootstrp(nBootstraps, @fitModel, data);

out=struct;
out.SpectralFit = normal_fit;
out.bootstat = bootstat;
out.bootsam = bootsam;

end

