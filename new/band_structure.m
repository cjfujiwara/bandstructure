
close all

npt=constants;
npt.depth=5;

[npt,hF1]=getBandStructure(npt);
[npt,hF2]=calculateGapTunneling(npt);
[npt,hF3]=wannier(npt);

hF4=showSwave;
hF5=showTrappedInteraction;
hF6 = show_lattice_shift_79;

hF7=showPwave;


%% Harmonic Ineraction
%{
Bvec=linspace(195,225,1000);

% Interaction energy parameters
g = @(as) 4*pi*hbar^2*as*a0/m;

% Trap angular frequency
trap_omega = @(U) 2*pi*sqrt(4*U)*fr;

% Harmoic oscillator length
trap_length = @(U) sqrt(hbar./(m*trap_omega(U)));

% Hubbard Interaction energy
U_hubbard_97 = @(U,B) g(feshbach_97(B))*(2*pi*trap_length(U).^2).^(-3/2);
U_hubbard_95 = @(U,B) g(feshbach_95(B))*(2*pi*trap_length(U).^2).^(-3/2);

str = ['$ U = \frac{4\pi\hbar^2 a}{m}\left(\frac{1}{2\pi \sigma_\mathrm{HO}^2}\right)^{3/2}$'];


hF6=figure(227);
set(hF6,'color','w','name','feshbach_energy');
hF6.Position=[420 540 700 400];
clf
co=get(gca,'colororder');
subplot(121);
p7_060=plot(Bvec,1e-3*U_hubbard_97(60,Bvec)/h,':','LineWidth',2,'color',co(1,:)); 
hold on
p7_100=plot(Bvec,1e-3*U_hubbard_97(100,Bvec)/h,'-.','LineWidth',2,'color',co(1,:)); 
p7_200=plot(Bvec,1e-3*U_hubbard_97(200,Bvec)/h,'-','LineWidth',2,'color',co(1,:)); 
xlim([195 202]);
ylim([0 300]);

xlabel('magnetic field (G)');
ylabel('interaction energy (kHz)');
set(gca,'box','on','xgrid','on','ygrid','on');

text(.98,.02,'$|ab\rangle$','interpreter','latex','horizontalalignment','right',...
    'verticalalignment','bottom','units','normalized','fontsize',14);


legend([p7_060,p7_100,p7_200],{'$60~E_R$','$100~E_R$','$200~E_R$'},'interpreter','latex',...
    'fontsize',10,'location','northwest');

text(.04,.78,str,'interpreter','latex','horizontalalignment','left',...
    'verticalalignment','top','units','normalized','fontsize',12);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');


subplot(122);

p5_060=plot(Bvec,1e-3*U_hubbard_95(60,Bvec)/h,'r:','LineWidth',2,'color',co(2,:)); 
hold on
p5_100=plot(Bvec,1e-3*U_hubbard_95(100,Bvec)/h,'r-.','LineWidth',2,'color',co(2,:)); 
p5_200=plot(Bvec,1e-3*U_hubbard_95(200,Bvec)/h,'r-','LineWidth',2,'color',co(2,:)); 
xlim([195 223]);

% ylim([0 300]);
set(gca,'box','on','xgrid','on','ygrid','on');


xlabel('magnetic field (G)');
ylabel('interaction energy (kHz)');


text(.98,.02,'$|ac\rangle$','interpreter','latex','horizontalalignment','right',...
    'verticalalignment','bottom','units','normalized','fontsize',14);


legend([p5_060,p5_100,p5_200],{'$60~E_R$','$100~E_R$','$200~E_R$'},'interpreter','latex',...
    'fontsize',10,'location','northwest');

text(.04,.78,str,'interpreter','latex','horizontalalignment','left',...
    'verticalalignment','top','units','normalized','fontsize',12);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');

hF7=figure(228);
set(hF7,'color','w','name','feshbach_differential');
hF7.Position=[1140 540 400 400];
clf

p60=plot(Bvec,1e-3*(U_hubbard_97(60,Bvec)-U_hubbard_95(60,Bvec))/h,':','LineWidth',2,'color','k'); 
hold on
p100=plot(Bvec,1e-3*(U_hubbard_97(100,Bvec)-U_hubbard_95(100,Bvec))/h,'--','LineWidth',2,'color','k'); 

p200=plot(Bvec,1e-3*(U_hubbard_97(200,Bvec)-U_hubbard_95(200,Bvec))/h,'-','LineWidth',2,'color','k'); 
xlim([199 202]);
ylim([0 300]);

xlabel('magnetic field (G)');
ylabel('differential interaction energy (kHz)');
set(gca,'box','on','xgrid','on','ygrid','on');

text(.98,.02,'$E_{ab}-E_{ac}$','interpreter','latex','horizontalalignment','right',...
    'verticalalignment','bottom','units','normalized','fontsize',14);


legend([p60,p100,p200],{'$60~E_R$','$100~E_R$','$200~E_R$'},'interpreter','latex',...
    'fontsize',10,'location','northwest');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');

hF7b=figure(231);
set(hF7b,'color','w','name','oscillator_length');
hF7b.Position=[1140 540 400 400];
clf

p60=plot(U,trap_length(U)/a0,'-','LineWidth',2,'color','k'); 

xlim([10 200]);
ylim([0 2e3]);

xlabel('lattice depth (E_R)');
ylabel('trap length (a_0)');
set(gca,'box','on','xgrid','on','ygrid','on');

% legend([p60,p100,p200],{'$60~E_R$','$100~E_R$','$200~E_R$'},'interpreter','latex',...
%     'fontsize',10,'location','northwest');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
str = ['$ U = \frac{4\pi\hbar^2 a}{m}\left(\frac{1}{2\pi \sigma_\mathrm{HO}^2}\right)^{3/2}$'];

sqrt(hbar./(m*trap_omega(U)));

str = ['$\sigma_{\mathrm{HO}}=\sqrt{\frac{\hbar}{m\omega}}=\sqrt{\frac{\hbar}{m 2\pi E_R/h \sqrt{4U_0}}}$'];
text(.04,.02,str,'interpreter','latex','horizontalalignment','left',...
    'verticalalignment','bottom','units','normalized','fontsize',12);
%% Fermi - Hubbard Parameters

hF8=figure(229);
set(hF8,'color','w','name','hubbard');
hF8.Position=[740 50 1000 300];
clf

t_hubbard = abs(BW(:,1))*fr/4;
U_hubbard = U_hubbard_97(U,0)'/h;


subplot(131);
pS=plot(U,U_hubbard,'LineWidth',2); 
xlabel('lattice depth (E_R)');
ylabel('interaction energy (Hz)');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
xlim([0 15]);
str = ['$ U = \frac{4\pi\hbar^2 a}{m}\left(\frac{1}{2\pi \sigma_\mathrm{HO}^2}\right)^{3/2}$'];
text(.02,.98,str,'interpreter','latex','units','normalized','fontsize',14,...
    'verticalalignment','top');

subplot(132);
pS=plot(U,t_hubbard,'LineWidth',2); 
xlabel('lattice depth (E_R)');
ylabel('tunneling rate (Hz)');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
xlim([0 15]);
str = ['$ t = \Delta_s/4$'];
text(.02,.98,str,'interpreter','latex','units','normalized','fontsize',14,...
    'verticalalignment','top');


subplot(133);
pS=plot(U,U_hubbard./t_hubbard,'k-','LineWidth',2); 
xlim([0 15]);
xlabel('lattice depth (E_R)');
ylabel('U/t');
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');



%% Harmonic Ineraction

funcs=trapped_interaction;
Eshift=funcs{2};

Bvec=linspace(198,225,1000);

% Interaction energy parameters
g = @(as) 4*pi*hbar^2*as*a0/m;

% Trap angular frequency
trap_omega = @(U) 2*pi*sqrt(4*U)*fr;

% Harmoic oscillator length
trap_length = @(U) sqrt(hbar./(m*trap_omega(U)));

% Hubbard Interaction energy
U_hubbard_97_sat = @(U,B) Eshift(a0*feshbach_97(B)./trap_length(U))*trap_omega(U)*hbar;
U_hubbard_95_sat = @(U,B) Eshift(a0*feshbach_95(B)./trap_length(U))*trap_omega(U)*hbar;


U_hubbard_97 = @(U,B) g(feshbach_97(B))*(2*pi*trap_length(U).^2).^(-3/2);
U_hubbard_95 = @(U,B) g(feshbach_95(B))*(2*pi*trap_length(U).^2).^(-3/2);

hFnew=figure(232);
clf
hFnew.Color='w';
hFnew.Position=[200 200 1200 400 ];
hFnew.Name='interaction_saturation_2';

depth=200;

subplot(131);
plot(Bvec,1e-3*U_hubbard_97(depth,Bvec)/h,'linewidth',2);
hold on
plot(Bvec,1e-3*U_hubbard_97_sat(depth,Bvec)/h,'linewidth',2);

plot([Bvec(1) Bvec(end)],1e-3*trap_omega(depth)/(2*pi)*[1 1],'k--','linewidth',1);

xlim([195 202]);
set(gca,'xgrid','on','ygrid','on','fontsize',12,'box','on','linewidth',1);
xlabel('magnetic field (G)');
ylim([0 400]);
xlim([198 202]);
ylabel('$U_{97} \mathrm{(kHz)}$','interpreter','latex');
legStr={['$' num2str(depth) '~E_R$, linear'],['$' num2str(depth) '~E_R$, saturated'],'band gap'};
legend(legStr,'interpreter','latex','location','northwest');


subplot(132);
plot(Bvec,1e-3*U_hubbard_95(depth,Bvec)/h,'linewidth',2);
hold on
plot(Bvec,1e-3*U_hubbard_95_sat(depth,Bvec)/h,'linewidth',2);
xlim([195 202]);
set(gca,'xgrid','on','ygrid','on','fontsize',12,'box','on','linewidth',1);
xlabel('magnetic field (G)');
xlim([198 202]);
ylabel('$U_{95} \mathrm{(kHz)}$','interpreter','latex');
legStr={['$' num2str(depth) '~E_R$, linear'],['$' num2str(depth) '~E_R$, saturated']};
legend(legStr,'interpreter','latex','location','northwest');

subplot(133);
plot(Bvec,1e-3*(U_hubbard_97(depth,Bvec)-U_hubbard_95(depth,Bvec))/h,'linewidth',2);
hold on
plot(Bvec,1e-3*(U_hubbard_97_sat(depth,Bvec)-U_hubbard_95_sat(depth,Bvec))/h,'linewidth',2);
plot([Bvec(1) Bvec(end)],1e-3*trap_omega(depth)/(2*pi)*[1 1],'k--','linewidth',1);
xlim([195 202]);
set(gca,'xgrid','on','ygrid','on','fontsize',12,'box','on','linewidth',1);
xlabel('magnetic field (G)');
ylim([0 150]);
xlim([198 202]);
ylabel('$\mathrm{U}_{97} - \mathrm{U}_{95} \mathrm{(kHz)}$','interpreter','latex');
legStr={['$' num2str(depth) '~E_R$, linear'],['$' num2str(depth) '~E_R$, saturated'],'band gap'};
legend(legStr,'interpreter','latex','location','southeast');
%}
%%
doSave=1;
if doSave
    disp('saving figures ...');    
    figs=get(groot,'children');
    
    for kk=1:length(figs)
        fprintf(['saving figure ' num2str(kk) ' ... ']);
        print(figs(kk),['figs/' figs(kk).Name],'-dpng','-r400'); 
        savefig(figs(kk),['figs/' figs(kk).Name]); 
        disp('done');

    end    
end
