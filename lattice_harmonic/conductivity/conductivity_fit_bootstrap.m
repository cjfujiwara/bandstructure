function out=conductivity_fit_bootstrap(freq,sigma,sigma_err)

freq=freq(:);
sigma=sigma(:);
sigma_err=sigma_err(:);

data=[freq sigma sigma_err];

%% Calculate Lattice Properties

lattice                     = constants;
lattice.depth               = [2.5]; 
lattice.WannierBands        = [1];

lattice.numStates           = 101;       % must be odd
lattice.numK                = 301;      % must be odd    

wannier_opts                = struct;
wannier_opts.Bands          = [1];


lattice = calculateBandStructure(lattice);      % calculate band structure
lattice = calculateTunneling(lattice);          % calculate tunneling elements
lattice = wannier(lattice,wannier_opts);        % Calculate wannier function
lattice = calculateWannierMoments(lattice);     % Dipole matrix elements in wannier basis

%% Fit it once

normal_fit = conductivity_fit2(lattice,freq,sigma,sigma_err);

%%
n=0;
     function fittedParams=fitModel(data)
        t1=now;
         freq = data(:,1);
        sigma = data(:,2);
        sigma_err = data(:,3);         
        output = conductivity_fit2(lattice,freq,sigma,sigma_err,normal_fit);
        fittedParams= [output.fout output.rho0 output.rhoinf]; 
        n=n+1;

        t2=now;
        fprintf([num2str(n) ' ' num2str(24*60*60*(t2-t1),'%.2f') ' sec.']);
        disp([num2str(fittedParams(1),'%.2f') ', ' num2str(fittedParams(2),'%.2f') ', ' num2str(fittedParams(3),'%.2f')]);

     end


nBootstraps = 1000;
% Apply bootstrap
[bootstat, bootsam] = bootstrp(nBootstraps, @fitModel, data);

out=struct;
out.SpectralFit = normal_fit;
out.bootstat = bootstat;
out.bootsam = bootsam;

end

