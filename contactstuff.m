bcs = @(x) 512/(945*pi^2).*x.^(-2).*(1+(256/(35*pi^2)-(63+189*log(2))/(1024))./x);
x1 = linspace(-200,-1,1e5);
ybcs = bcs(x1);

bec = @(x) x + 5^(2/5)/(2^(12/5)*7)*0.6^(2/5).*x.^(-7/5);
x2 = linspace(0.5,30,1e5);
ybec = bec(x2);


U = 0.28;

xq = linspace(-2,1,1e3);
s = spline([x1 0 x2],[ybcs U ybec],xq);

Finterp = @(xq) spline([x1 0 x2],[ybcs U ybec],xq);


%% Plot it
hF=figure(10);
hF.Position=[100 100 400 350];
clf
set(gcf,'color','w');
co=get(gca,'colororder');
hold on
pBCS=plot(x1,ybcs,'-','linewidth',2,'color',co(1,:));
hold on
pBEC=plot(x2,ybec,'-','linewidth',2,'color',co(2,:));

pS=plot(xq,s,'--','linewidth',1,'color','k');

pU=plot(0,U,'o','linewidth',2,'color','k','markerfacecolor','k');
xlabel('$1/(k_F a)$','interpreter','latex');
ylabel('$\mathcal{F}$','interpreter','latex');
legend([pBCS pBEC pU pS],{'BCS','BEC','unitary','spline'},'location','northwest');
set(gca,'fontsize',14,'fontname','times','box','on',...
    'xgrid','on','ygrid','on','linewidth',1);
ylim([-.05 1.6]);
xlim([-1.5 1.5]);


hF=figure(11);
hF.Position=[100 100 400 350];
clf
set(gcf,'color','w');
co=get(gca,'colororder');
hold on
pBCS=plot(x1,ybcs,'-','linewidth',2,'color',co(1,:));
hold on
pBEC=plot(x2,ybec,'-','linewidth',2,'color',co(2,:));

pS=plot(xq,s,'--','linewidth',1,'color','k');

pU=plot(0,U,'o','linewidth',2,'color','k','markerfacecolor','k');
xlabel('$1/(k_F a)$','interpreter','latex');
ylabel('$\mathcal{F}$','interpreter','latex');
legend([pBCS pBEC pU pS],{'BCS','BEC','unitary','spline'},'location','northwest');
set(gca,'fontsize',14,'fontname','times','box','on',...
    'xgrid','on','ygrid','on','linewidth',1);
ylim([.5*1e-2 2]);
xlim([-3 2]);
set(gca,'YScale','log');


%% 

% constants
amu = 1.66054e-27; % amu in kg
h = 6.62607015e-34; % planck's constant  in Js
hbar = h/(2*pi); % reduced planck's constant Js
a0 = 5.29177e-11;       % bohr radisu in m
mubh = 1.39962449e6; % bohr magneton/h in Hz/Gauss
mub = mubh*h; % (Bohr magneton in J/Gauss).
kB = 1.381e-23 ; % boltzmann constant in J/K

m = 40*amu; % amss
mu = 1.5*mub; % magnetic fmoment

% Fermig as
fbar = 100; % trap frequency in Hz
N=1e5;
Ef = (3*N)^(1/3)*(h*fbar);
Tf = Ef/kB;
kF = sqrt(2*m*Ef/hbar^2);

% atom


% feshbach
a_bg = 166.978*a0;
Delta = 6.910;
B0 = 202.15;

% Feshbach field
B2a = @(B) a_bg*(1-Delta./(B-B0));


s1 = ['$N_b = Nk_F R_* \mathcal{F}(1/(k_F a))(1-a_{bg}/a)^2$' newline ...
    '$R = \frac{\hbar^2}{m a_\mathrm{bg}\mu  \Delta}$'];

s2 = ['$N = 10^5, T_F = ' num2str(round(Tf*1e9,1)) '~\mathrm{nK}$'];

% Moleucloar raidus thingamajigger
Rstar = hbar^2./(m*a_bg*mu*Delta);
Bvec = linspace(190,208,1e3);
Nb = N*kF*Rstar*Finterp(1./(kF*B2a(Bvec))).*(1-a_bg./B2a(Bvec)).^2;

%%
hF2 = figure(12)
hF2.Color='w';

hF2.Position = [100 100 400 400];
clf
plot(Bvec,Nb,'k-','linewidth',2);
xlabel('field (G)');
ylabel('N_b');
text(.02,.02,s1,'interpreter','latex','verticalalignment','bottom',...
    'fontsize',12,'units','normalized');
text(.98,.98,s2,'interpreter','latex','verticalalignment','top',...
    'fontsize',12,'units','normalized','horizontalalignment','right');
set(gca,'fontsize',14,'fontname','times','box','on','linewidth',1,...
    'xgrid','on','ygrid','on');
xlim([202 208]);

% set(gca,'YScale','log');
yL = get(gca,'Ylim');
ylim([0 yL(2)]);

 yyaxis right
plot(Bvec,(1./(kF*B2a(Bvec))),'linewidth',2)
ylabel('1/(k_Fa)')

%%

Omega = 2*pi*1.3e6;
gamma = 2*pi*26e6;


R0 = 2*(Nb/N)*Omega^2/gamma;

s2 = ['$N = 10^5, T_F = ' num2str(round(Tf*1e9,1)) '~\mathrm{nK}$' newline ...
    '$\Omega = 2\pi\cdot' num2str(1e-6*Omega/(2*pi)) '~\mathrm{MHz}$, ' ...
    '$\gamma = 2\pi\cdot' num2str(1e-6*gamma/(2*pi)) '~\mathrm{MHz}$'];

hF2 = figure(13);
hF2.Color='w';

hF2.Position = [100 100 400 400];
clf
plot(Bvec,R0*1e-3,'k-','linewidth',2);
xlabel('field (G)');
ylabel('Loss Rate (kHz)');
text(.02,.02,s1,'interpreter','latex','verticalalignment','bottom',...
    'fontsize',12,'units','normalized');
text(.98,.98,s2,'interpreter','latex','verticalalignment','top',...
    'fontsize',12,'units','normalized','horizontalalignment','right');
set(gca,'fontsize',14,'fontname','times','box','on','linewidth',1,...
    'xgrid','on','ygrid','on');
xlim([202 208]);

% set(gca,'YScale','log');
yL = get(gca,'Ylim');
ylim([0 yL(2)]);
% 
%  yyaxis right
% plot(Bvec,(1./(kF*B2a(Bvec))),'linewidth',2)
% ylabel('1/(k_Fa)')
