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
nBootstraps = 1000;


D = parallel.pool.DataQueue;
h = waitbar(0,'Please wait ...');
afterEach(D,@nUpdateWaitbar);

N = 200;
p = 1;


n=0;

     function fittedParams=fitModel(data)
         D = parallel.pool.DataQueue;        

        t1=now;
         freq = data(:,1);
        sigma = data(:,2);
        sigma_err = data(:,3);         
        output = conductivity_fit2(lattice,freq,sigma,sigma_err,normal_fit);
        fittedParams= [output.fout output.rho0 output.rhoinf]; 

        t2=now;
        % fprintf([num2str(n) ' ' num2str(24*60*60*(t2-t1),'%.2f') ' sec.']);
        % fprintf([num2str(24*60*60*(t2-t1),'%.2f') ' sec ']);
        disp([num2str(fittedParams(1),'%.2f') ', ' num2str(fittedParams(2),'%.2f') ', ' num2str(fittedParams(3),'%.2f')]);
        % n=n+1;
        send(D,0);
     end

 function nUpdateWaitbar(~)
        waitbar(p/nBootstraps,h);
        p = p + 1;
    end



options.UseParallel	=true;
options.UseSubstreams	=false;
% options.UseParallel	='true';



% Apply bootstrap
tic
[bootstat, bootsam] = bootstrp(nBootstraps, @fitModel, data,'Options',options);
toc
out=struct;
out.SpectralFit = normal_fit;
out.bootstat = bootstat;
out.bootsam = bootsam;

end

