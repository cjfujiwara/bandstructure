%% Introduction
% This script calculates the properties of an optical lattice with an
% additional harmonic confinement.
%


%% Initialize
% Define parameters of calculation

lattice=constants;
lattice.depth=[2.5]; 

%% Flags
doShowBandStructure = true;
doAnimateWannier = false;

%% Caclulate Band Properties
% Calculate the band structure

lattice = calculateBandStructure(lattice);   % calculate band structure

if doShowBandStructure
    % Plot the band structure
    show_band_opts = struct;
    show_band_opts.Bands = 1:3;
    hF_band = showBandStructure(lattice,show_band_opts);
end

%% Calculate Tunneling Propertiess
% Calculate the tunneling matrix elements
lattice = calculateTunneling(lattice);      % calculate tunneling elements

%% Calculate Wannier
% Calculate the wannier functions, specify which bands you want to
% calculate
wannier_opts = struct;
wannier_opts.Bands = [1];
lattice.WannierBands = wannier_opts.Bands;

lattice = wannier(lattice,wannier_opts);                % Calculate wannier function
lattice = calculateWannierMoments(lattice);             % Dipole matrix elements in wannier basis

% Show the Wannier function
hF_wannier = showWannier(lattice,wannier_opts);           % calculate wannier function 
%% Wannier Harmonic Coupling
% Calculate the matrix coupling element induced from a harmonic potential
% on the wannier states (this is primarily important for multi-band
% physics)
%
% <w_m(x_i)|x^2|w_n(x_j)>
% <w_m(x_i)|x^1|w_n(x_j)>

lattice=calculateWannierHarmonicCoupling2(lattice);
%% Wannier Animation
% Animate the wannier functions if you specified different lattice depths
if doAnimateWannier
    tempfile = fullfile(tempdir,'animate.gif');
    for kk=1:length(hF_wannier)    
        frame = getframe(hF_wannier(kk));
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);  
        if kk == 1
            imwrite(A,map,tempfile,'gif','LoopCount',Inf,'DelayTime',1);
        else
            if kk==length(hF_wannier)
                imwrite(A,map,tempfile,'gif','WriteMode','append','DelayTime',1);
            else
                imwrite(A,map,tempfile,'gif','WriteMode','append','DelayTime',.1);
            end
        end        
    end
    copyfile(tempfile,'wannier.gif','f');
end

%% Harmonic Coupling

 lattice = calculateWannierHarmonicCoupling(lattice);

%% Calculate 1D spectrum with Harmonic Confinement

% calculation parameters
harmonic_opts = struct;
harmonic_opts.NumSites =601;
harmonic_opts.MaxTunnelingOrder = 51;
harmonic_opts.NumBands =1;



% XY Lattice
harmonic_opts.omega = 2*pi*57;
harmonic_opts.omega = 2*pi*55;

harmonic_opts.Omega = 0.5*lattice.m*harmonic_opts.omega^2*(lattice.lambda/2)^2/lattice.h;
[lattice,harmonic_output_H] = calculateLHOSpectrum_sband(lattice,harmonic_opts);
% [lattice,harmonic_output_H2] = calculateLHOSpectrum(lattice,harmonic_opts);

% Z Direction
harmonic_opts.omega = 2*pi*266; % XDT Vertical trap frequency
harmonic_opts.Omega = 0.5*lattice.m*harmonic_opts.omega^2*(lattice.lambda/2)^2/lattice.h;
% [lattice,harmonic_output_V] = calculateLHOSpectrum(lattice,harmonic_opts);
[lattice,harmonic_output_V] = calculateLHOSpectrum_sband(lattice,harmonic_opts);

% Fit lowest band to linear dispersion
harmonic_output_H = fitHOtoFirstBand(harmonic_output_H);
harmonic_output_V = fitHOtoFirstBand(harmonic_output_V);

%%

hF_x=showLatticeHarmonic(harmonic_output_H,lattice);
xlim([0 60]);
ylim(-6500 + [0 3000])

hF_z=showLatticeHarmonic(harmonic_output_V,lattice);
xlim([0 20]);
ylim(-6500 + [0 3000])
hF_z.Position(1) = hF_x.Position(1)+hF_x.Position(3)+5;

%% Eigen States
showLHO_Eigenstates(lattice,harmonic_output_H)
showLHO_Eigenstates(lattice,harmonic_output_V)



%% Show Differential Energy
out=harmonic_output_H;
hF_eng_diff = figure(1010);
clf
hF_eng_diff.Color='w';

uu=1;

ax1 = axes;
% Actual Energy Data
[~,dominateBandIndex] = max(out.BandProjection(:,:,uu),[],2);
co=get(gca,'colororder');
colors = co(mod(dominateBandIndex-1,7)+1,:);
inds = 1:size(out.EigenValues(:,uu));
E_min = min(out.EigenValues(:,uu));
nstates = out.NumSites*out.NumBands;
eng = out.EigenValues(:,uu);

E0=eng(1);
pData=scatter(eng(1:end-1)-E0,diff(eng),2,colors(1:end-1,:),'linewidth',2,...
    'parent',ax1);
ylim([0 130])
xlim([0 20]*563);
ylabel('${E}_{n+1}-{E}_{n}$ [Hz]','interpreter','latex')
hold on

t=563;
T=[1 2 3 4 5 6 7 8]*t;
f_vec=linspace(0,30*t,100);

myc = jet(length(T));
set(gca,'box','on','linewidth',1,'fontsize',12)

xlabel('$E_n-E_0$ [Hz]','interpreter','latex')
ax2 = axes;
ax2.Position=ax1.Position;

clear ps
for nn=1:length(T)
    ps(nn)=plot(f_vec,exp(-f_vec/T(nn)),'-','color',[myc(nn,:) .5]);
    strs{nn}=['T/t = ' num2str(T(nn)/t)];
    hold on
end
ylim([0 1])
xlim([0 20]*t)
set(ax2,'Visible','off')
linkaxes([ax1 ax2],'x');

legend([ps(1) ps(end)],strs([1 length(T)]),'location','southeast');

%%
% showBandProjections(harmonic_output_H)
% showBandProjections(harmonic_output_V)

%% 
opts=struct;
opts.Indeces = 'auto';
opts.Indeces = [1 2 50 51];
% opts.Indeces = [1:100];

showLatticeHarmonicWavefunction(lattice,harmonic_output_H,opts);
%%
opts.Indeces = [1:150];
wfs = calculateLatticeHarmonicWavefunction(lattice,harmonic_output_H,opts);
%%
D=zeros(size(wfs,2),size(wfs,2));
for rr=1:size(wfs,2)
    for cc = 1:size(wfs,2)
        D(rr,cc)=trapz(conj(wfs(:,rr)).*wfs(:,cc).*x);
    end
end
%%


N=size(wfs,2);

out=harmonic_output_H;
[~,dominateBandIndex] = max(out.BandProjection(:,:,uu),[],2);
co=get(gca,'colororder');
colors = co(mod(dominateBandIndex-1,7)+1,:);
inds = 1:size(out.EigenValues(:,uu));
E_min = min(out.EigenValues(:,uu));
nstates = out.NumSites*out.NumBands;
eng = out.EigenValues(:,uu);


figure(999);

subplot(1,2,1)
imagesc(abs(real(D)));set(gca,'YDir','normal');colorbar
set(gca,'fontsize',10);
xlabel('eigen index');
ylabel('eigen index');
axis equal tight
title('$|\langle \psi_m|x|\psi_n\rangle|$','interpreter','latex','fontsize',18)

subplot(1,2,2)
% plot(harmonic_output_H.EigenValues(1:200))

pData=scatter(inds(1:N),eng(1:N)-eng(1),2,colors(1:N,:),'linewidth',2);
set(gca,'box','on','linewidth',1,'fontsize',12)
xlabel('eigenindex');
ylabel('energy - E_0 [Hz]');
title('eigenspectrum 2.5 Er + 60 Hz HO')
%% Thermodynamical Analysis

% calculateThermodynamics(lattice,harmonic_output_H,harmonic_output_H,harmonic_output_V);



 