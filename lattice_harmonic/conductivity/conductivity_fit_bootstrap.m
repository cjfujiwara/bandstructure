function out=conductivity_fit_bootstrap(bs_moment)


freq            = [bs_moment.Frequency_Hz];
sigma           = [bs_moment.sigma];
sigma_err       = [bs_moment.sigmaErr];

% Interpreter R2 as an "error" for weighting purposes
R2      = [bs_moment.rsquare];
R2_err  = sqrt(1./R2);
R2_err  = (1+1i)*R2_err.*abs(mean(real(sigma_err)));        
sigma_err = R2_err;


% out(nn)=conductivity_fit_bootstrap(f,s,s_err);

varX_um = arrayfun( @(i) mean(bs_moment.VarianceX_um{i}),1:length(bs_moment.VarianceX_um));
varY_um = arrayfun( @(i) mean(bs_moment.VarianceY_um{i}),1:length(bs_moment.VarianceY_um));


varX_um2 = varX_um(:);
varY_um2 = varY_um(:);

freq        = freq(:);
sigma       = sigma(:);
sigma_err   = sigma_err(:);
data        = [freq sigma sigma_err varX_um2 varY_um2];

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
nBootstraps =100;
n=0;

     function fittedParams=fitModel(data)
        t1=now;
        F = data(:,1);
        S = data(:,2);
        SE = data(:,3);  
        X2 = data(:,4);
        Y2 = data(:,5);
        output = conductivity_fit2(lattice,F,S,SE,normal_fit);
        fittedParams= [output.fout output.rho0 output.rhoinf]; 
        t2=now;
      
        amu= 1.66054e-27; %[kg];
        % m = 40 * amu * lattice.BandMassGamma(1,1);
        m = 40 * amu;

        h = lattice.h;

        Tx = m*mean(X2)*(1e-12)*(2*pi*output.fout(3))^2/h;
        Ty = m*mean(Y2)*(1e-12)*(2*pi*output.fout(3))^2/h;

% keyboard
        fittedParams = [fittedParams Tx Ty];

          fprintf([num2str(24*60*60*(t2-t1),'%.2f') ' sec ']);
            disp([num2str(fittedParams(1),'%.2f') ', ' ...
                num2str(fittedParams(2),'%.2f') ', ' ...
                num2str(fittedParams(3),'%.2f') ', ' ...
                num2str(fittedParams(6),'%.2f') ', ' ...
                num2str(fittedParams(7),'%.2f') ]);
            n=n+1;
     end

options.UseParallel	=true;
options.UseSubstreams	=false;

    FREQ_THEORY = linspace(0,200,200);


    figure(20);
    clf
    co=get(gca,'colororder');
    errorbar(freq,real(sigma),real(sigma_err),'o','markerfacecolor',co(1,:));
    hold on
    errorbar(freq,imag(sigma),imag(sigma_err),'o','markerfacecolor',co(2,:));
    drawnow;
    yF=conductivity_eval2(FREQ_THEORY, [normal_fit.fout],lattice);
    plot(FREQ_THEORY,real(yF),'-','color',co(1,:));
    hold on
    plot(FREQ_THEORY,imag(yF),'-','color',co(2,:));
    drawnow;

% keyboard

% Apply bootstrap
tic
[bootstat, bootsam] = bootstrp(nBootstraps, @fitModel, data,'Options',options);
toc
out=struct;
out.SpectralFit = normal_fit;
out.bootstat = bootstat;
out.bootsam = bootsam;

end

