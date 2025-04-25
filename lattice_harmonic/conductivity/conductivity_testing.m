T=1126;
G=2*pi*20;

ff = 5:1:400;
zz=conductivity2(ff,T,G,lho);
rho = 1./zz;

u = lho.Depth;
t = [lho.Tunneling];
Nt = lho.MaxTunnelingOrder;
f = lho.omega/(2*pi);

str=['$V_0=' num2str(u) 'E_r$' newline ...
    '$t=[' num2str(round(t(1))) ...
    ',' num2str(round(t(2))) ',\dots]~\mathrm{Hz},(\Delta j_\mathrm{max}=' num2str(Nt) ')$'];

str1 = ['$(T=' num2str(round(T/t(1),1)) 't,' ...
    '\Gamma=' num2str(G) '/\mathrm{s},' ...
    'f_\mathrm{HO}=' num2str(f) '~\mathrm{Hz})$' ];

hF=figure(5);
clf
hF.Color='w';
hF.Position=[50 50 800 350];

subplot(121);
plot(ff,real(zz))
hold on
plot(ff,imag(zz))
xlabel('frequency (Hz)');
ylabel('conductivity (\sigma_0)');
legend({'Re','Im'});
text(.99,.01,str,'interpreter','latex','units','normalized',...
    'horizontalalignment','right','verticalalignment','bottom');
title(['conductivity ' str1],'interpreter','latex');

subplot(122);
plot(ff,real(rho))
hold on
plot(ff,imag(rho))
xlabel('frequency (Hz)');
ylabel('resitivity (1/\sigma_0)');
legend({'Re','Im'});
title(['resitivity ' str1],'interpreter','latex');
text(.99,.01,str,'interpreter','latex','units','normalized',...
    'horizontalalignment','right','verticalalignment','bottom');
ylim([-.1 .1]);
%%

data    = load('C:\Users\coraf\Downloads\2025.03.16 198.5 G spectrum data.mat');
data    = data.Mar16spectrum198p5;
z       = data.sR+ 1i*data.sI;
f       = data.freq;
fout    = conductivity_fit(f,z);

%%
inds=[1:15];


clear myfits
for kk=1:length(inds)
    x = composite_data.digdata(kk).X;
    y = composite_data.digdata(kk).Xc_um;
    
    P = [composite_data.digdata(kk).Params];
    T = [P.conductivity_mod_time];
    Tr = [P.conductivity_mod_ramp_time];
    Ttot = T+Tr;
    
    
    fme = [composite_data.digdata(kk).Params.conductivity_mod_freq];
    fme = fme(1)*1e-3;
    
    myfits(kk)=bootstrap_oscillations_fit(Ttot,y,fme);
   
end
