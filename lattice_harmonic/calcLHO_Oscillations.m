function calcLHO_Oscillations


lattice                     = constants;
lattice.depth               = [2.5]; 
lattice.WannierBands        = [1];

lattice.numStates           = 101;       % must be odd
lattice.numK                = 301;      % must be odd    

wannier_opts                = struct;
wannier_opts.Bands          = [1];

Nsites = 201;
TunnelOrder = 11;
HarmonicBands =1 ;

lattice = calculateBandStructure(lattice);      % calculate band structure
lattice = calculateTunneling(lattice);          % calculate tunneling elements
lattice = wannier(lattice,wannier_opts);        % Calculate wannier function
lattice = calculateWannierMoments(lattice);     % Dipole matrix elements in wannier basis


m0star = lattice.BandMassGamma(1); % band mass in units of bare mass

displacement_sites = 1;     % 10 sites is like 5 um
TrapFrequency_Hz = 45:.5:80;
Temperature_Hz   = 0:50:2500;
Temperature_Hz(1) =0.01;

[fff,ttt]=meshgrid(TrapFrequency_Hz,Temperature_Hz);

gauss_radii      = zeros(size(fff,1),size(fff,2));
osc_freq        = zeros(size(fff,1),size(fff,2));
osc_tau         = zeros(size(fff,1),size(fff,2));
osc_cen         = zeros(size(fff,1),size(fff,2));
osc_amp         = zeros(size(fff,1),size(fff,2));


dt = 1e-3;
tVec = 0:dt:0.2;

freq_fit_list=30:1:100;
corr_mat=zeros(length(freq_fit_list),length(tVec));
for a =1:length(freq_fit_list)
    corr_mat(a,:) = cos(2*pi*freq_fit_list(a)*tVec);
end




%%
decay_cosine = fittype('A*cos(2*pi*f*t)*exp(-t/tau)+B','independent','t',...
    'coefficients',{'A','tau','f','B'});
fit_opts = fitoptions(decay_cosine);
%%


for ff=1:length(TrapFrequency_Hz)
    fprintf([num2str(ff) '/' num2str(length(TrapFrequency_Hz)) ' ']);
    fprintf('solving eigenvalue...')
    % Solve Eigenvalue problem
    omega = 2*pi*TrapFrequency_Hz(ff);
    opts=struct;               
    opts.NumSites = Nsites;
    opts.MaxTunnelingOrder = TunnelOrder;
    opts.HarmonicBands = HarmonicBands;
    opts.omega = omega;
    opts.Omega = 0.5*lattice.m*omega^2*(lattice.lambda/2)^2/lattice.h;
    [~,lho] = calculateLHOSpectrum_sband(lattice,opts);

    % Formatting
    probEig = abs(lho.EigenVectors).^2;     % probabiity amplitude of each eigenstates
    eng_Hz = lho.EigenValues-lho.EigenValues(1); % recast energy
    X= lho.PositionVector;
    X = X(:);
    H=lho.Hamiltonian-eye(Nsites)*lho.EigenValues(1);

    % Time Evolution Operator
    U_exp = expm(-1i*H*dt*2*pi);% exponentiaion time evolution
    % U_cn = (1-1i*0.5*H*dt)/(1+1i*0.5*H*dt); % implicit crank nicolson

    CoM_all = zeros(length(Nsites),length(tVec));     
    fprintf(' oscillating ...')
    for pp=1:Nsites        
        % fprintf([num2str(pp) ','])
        psi = lho.EigenVectors(:,pp);
        psi = circshift(psi,displacement_sites);
        % [t,y] = ode45(@(t,y) 1i*H*y,0:1e-3:0.2,psi,options);
        y = zeros(length(psi),length(tVec));
        y(:,1)=psi;
        for n=2:length(tVec)
            y(:,n)= U_exp*y(:,n-1);
        end
        Xall=repmat(X,[1 size(y,2)]);
        CoM = sum(Xall.*abs(y).^2,1);
        CoM_all(pp,:)=CoM;            
    end

     fprintf(' temperature ...')

    figure(6);
    clf
    for tt=1:length(Temperature_Hz)
        Z = sum(exp(-eng_Hz/Temperature_Hz(tt)));  % Partition Function
        prob = exp(-eng_Hz/Temperature_Hz(tt))/Z;  % Eigenstate probability
        prob = prob';                              % Transpose
        probMat = repmat(prob,[Nsites 1]);
        probTot = sum(probMat.*probEig,2);         % prbability density
        sigma = sqrt(sum(probTot.*X.^2));          % Second moment
        gauss_radii(tt,ff)=sigma;
        CoM_temp = CoM_all.*repmat(prob(:),[1 size(CoM_all,2)]);
        CoM_temp = sum(CoM_temp,1);   

        [~,ind]=max(sum(corr_mat.*repmat(CoM_temp,[size(corr_mat,1) 1]),2));
        f_guess = freq_fit_list(ind);
        figure(6);
        plot(tVec,CoM_temp)
        hold on
        fit_opts.StartPoint = [displacement_sites tVec(end) f_guess mean(CoM_temp)];
        fit_opts.Upper = [10 inf 200 10];
        fit_opts.Lower = [0 .01 30 -1];

        Tmax = 3/f_guess;

        W=[tVec(:)<=Tmax];
        fit_opts.Weights = double(W);
        
        fout=fit(tVec(:),CoM_temp(:),decay_cosine,fit_opts);
        osc_freq(tt,ff) = fout.f;
        osc_tau(tt,ff) = fout.tau;
        osc_cen(tt,ff) = fout.B;
        osc_amp(tt,ff) = fout.A;     
    end      
    disp(' done')
end

%% Show Results
hF_size = figure;
hF_size.Color='w';
clf

subplot(141);
imagesc(TrapFrequency_Hz,Temperature_Hz,gauss_radii*0.532);
set(gca,'YDir','normal','fontsize',14,'fontname','times');
xlabel('trap frequency (Hz)');
ylabel('temperature (Hz)')
cc=colorbar;
cc.Label.String = 'second moment (um)';
title('size')

subplot(142);
imagesc(TrapFrequency_Hz,Temperature_Hz,osc_freq);
set(gca,'YDir','normal','fontsize',14,'fontname','times');
xlabel('trap frequency (Hz)');
ylabel('temperature (Hz)')
cc=colorbar;
cc.Label.String = 'oscillation frequency (Hz)';
title('frequency')



m0star = lattice.BandMassGamma(1); % band mass in units of bare mass
subplot(144);
imagesc(TrapFrequency_Hz,Temperature_Hz,sqrt(m0star)*osc_freq./fff);
set(gca,'YDir','normal','fontsize',14,'fontname','times');
xlabel('trap frequency (Hz)');
ylabel('temperature (Hz)')
cc=colorbar;
cc.Label.String = '$\sqrt{m_0^*/m_0}f_\mathrm{osc}/f_\mathrm{trap}$';
cc.Label.Interpreter='latex';
title('normalized oscillation frequency')
caxis([.9 1]);


subplot(143);
cla
myslist = [6:.5:8];

for ss = 1:length(myslist)
    s=myslist(ss);
    T0=zeros(length(TrapFrequency_Hz),1);
    freq=T0;
    for ff = 1:length(TrapFrequency_Hz)
        this_s_list = gauss_radii(:,ff)*0.532;
        this_freq = osc_freq(:,ff);
        T0(ff) = fzero(@(T) interp1(Temperature_Hz,this_s_list,T)-s ,1000);
        freq(ff) = interp1(Temperature_Hz,this_freq,T0(ff));
    end
    plot(freq,TrapFrequency_Hz(:)./freq);
    hold on
    legStr{ss}=['$\sigma =' num2str(s) '\mu\mathrm{m}$'];
end
set(gca,'YDir','normal','fontsize',14,'fontname','times');
legend(legStr,'interpreter','latex','fontsize',10,'location','southeast');

xlabel('oscillation frequency (Hz)');
ylabel('trap frequency/oscillation frequency')
% keyboard



end

function [t, y] = RK4(H,psi0,tf,dt)
    % Initialize time vector
    N = round(tf/dt);
    t = linspace(0,tf,N);
    y = zeros(numel(psi0),N);
    y(:,1)= psi0;

    % Runge-Kutta 4th order loop
    for ii = 1:(N-1)
        k1 = 1i*H*y(:,ii);
        k2 = 1i*H*(y(:,ii)+0.5*dt*k1);
        k3 = 1i*H*(y(:,ii)+0.5*dt*k2);
        k4 = 1i*H*(y(:,ii)+dt*k3);        
        y(:,ii+1) = y(:,ii) + (k1 + 2*k2 + 2*k3 + k4) / 6;
        
    end
    % keyboard
end