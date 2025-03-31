tic
close all
opts = struct;
opts.doPlot = 0;

npt=constants;
npt.depth=.1:.1:10;

[npt,hF1]=getBandStructure(npt,opts);
% [npt,hF1_a]=calculateAMCoupling(npt,opts);

%% Wannier Calculation
opts.doPlot=0;
 [npt,hF3]=wannier(npt,opts);
 opts.doPlot = 1;
 
 %% Overlap Integrals
 npt = wannier_overlap_s(npt);
 
 [npt,hF2]=calculateGapTunneling(npt);

 hFa = showSwave(npt);
 %%
 b = [npt.OverlapS];
 Uho = [b.U_Harmonic_Theory_Er];
  Uw = [b.U_Wannier_Er];
  Er = 4.48;
  u0=[npt.OverlapS.depth];
  
  
  ut = [npt.bandwidth.depth];
  t = [npt.bandwidth.bandwidth];
  %%
  f=figure(20);
  f.Color='w';
  clf
  
%   U := g Integral |w(x)|^4 d^3x = g \prod_i |w(x_i)|^4 dx_i
%
% where g = 4*pi*a*hbar^2/m
  
subplot(131);
  str =['$U = \frac{4\pi a\hbar^2}{m} \int \mathrm{d}x |w(x)|^4 $'] ;
  
  a = 167*npt.a0;
  a = a/(532e-9);
  aho = [npt.Harmonic_Length];  
  plot(u0,a./aho.*Uw(1,:)*Er,'linewidth',2);hold on;
    plot(u0,a./aho.*Uho(1,:)*Er,'linewidth',2);
    legend({'wannier','harmonic'},'location','southeast');
    xlabel('lattice depth');
    ylabel('interaction energy (kHz)');
    text(.05,.95,str,'units','normalized','interpreter','latex','fontsize',16,...
        'horizontalalignment','left','verticalalignment','top');
    set(gca,'fontsize',14,'xgrid','on','ygrid','on');
    
    subplot(132);
    plot(ut,t(:,1)*Er/4,'linewidth',2)
    xlim([0 10]);
    xlabel('lattice depth (Er)');
    ylabel('tunneling (kHz)');
        set(gca,'fontsize',14,'xgrid','on','ygrid','on');
        legend({'\Delta_{1D}/4','wanier overlap'});


    