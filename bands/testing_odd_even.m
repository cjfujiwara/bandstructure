% calculation parameters
harmonic_opts = struct;
harmonic_opts.NumSites =601;
harmonic_opts.MaxTunnelingOrder = 5;
harmonic_opts.NumBands =1;

% XY Lattice
harmonic_opts.omega = 2*pi*57;
harmonic_opts.Omega = 0.5*lattice.m*harmonic_opts.omega^2*(lattice.lambda/2)^2/lattice.h;

harmonic_opts.NumSites =301;
[~,bb]=calculateLatticeHarmonicSpectrumOdd(npt,harmonic_opts);
xp = 0:300;

harmonic_opts.NumSites =601;
[~,aa]=calculateLatticeHarmonicSpectrum(npt,harmonic_opts);
x = -300:300;
nc = 301;
%%

ind_p = 40;
figure(20);
clf

yp = bb.EigenVectors(:,ind_p)/sqrt(2);
y  = aa.EigenVectors(:,2*ind_p);

sp = sign(yp(1));
s = sign(y(nc+1));

if s==sp
   p=1;
else
    p=-1;
end

plot(xp,yp,'-','linewidth',2);
hold on
plot(x,y*p,'-','linewidth',2);

xlim([0 100])

%%

figure(21)
clf

ea = aa.EigenValues(2:2:100)-aa.EigenValues(2);
eb=bb.EigenValues(1:50)-bb.EigenValues(1);

subplot(121);
 plot(ea,'o'); hold on
plot(eb,'o')

subplot(122);
plot(ea-eb);



