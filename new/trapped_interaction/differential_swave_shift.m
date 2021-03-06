function differential_swave_shift

out     = constants;
a0      = out.a0;
fr      = out.fr;
U       = 200;

% Band gap
Delta   = sqrt(4*U)*fr;
% Delta   = [122529.554801461];

% Harmonic oscillator length scale
aho=harmonic_length(U)/a0;


Bvec = linspace(190,225,1e4);
aVec = linspace(-10,10,1e4);

% Calculat the resonance shift
[dB,~]=lattice_shift_79(U,U,U);

% Calculate the theoretical function
[funcs] = trapped_interaction;
int_func=funcs{2};

aLinVec = linspace(-.05,.05,100);
myfit=fittype('A*x','independent','x','coefficient','A');
opt=fitoptions(myfit);
opt.StartPoint = sqrt(2/pi);
fout1 = fit(aLinVec',int_func(aLinVec'),myfit,opt);
mLin = fout1.A;


% mLin = int_func(.05)/.05;

% keyboard

% Get the feshbach resonance
a97 = feshbach_97(Bvec+dB);
a95 = feshbach_95(Bvec+dB);
 
U_97 = Delta*int_func(a97/aho);
U_95 = Delta*int_func(a95/aho);

mLin = fout1.A;
mLin = sqrt(2/pi);

% U_97_linear = Delta*(a97/aho)*sqrt(2/pi);
% U_95_linear = Delta*(a95/aho)*sqrt(2/pi);

U_97_linear = Delta*(a97/aho)*mLin;
U_95_linear = Delta*(a95/aho)*mLin;



%% THeoery
%{
hF=figure(1);
clf

set(gcf,'color','w');
hF.Position=[50 50 800 800];

subplot(221);
% pFree=plot(aVec,aVec*sqrt(2/pi),'k--','linewidth',2);
pFree=plot(aVec,aVec*mLin,'k--','linewidth',2);


hold on
co=get(gca,'colororder');
for kk=1:4
    ps(kk)=plot(aVec,funcs{kk}(aVec),'-','linewidth',2,...
        'color',co(kk,:));
    hold on
end



% legStr={'bound','first unbound','linear'};
% legend([ps(1:2) pFree],legStr,'location','southwest');

ylim([-3 3]);
xlim([-2 2]);

set(gca,'xgrid','on','ygrid','on','fontsize',10,'box','on',...
    'linewidth',1);

xlabel('scattering length $(a_\mathrm{HO})$','interpreter','latex');
ylabel('binding energy ($\hbar \omega$)','interpreter','latex');

subplot(222);
% pFree=plot(aVec,aVec*sqrt(2/pi),'k--','linewidth',2);
pFree=plot(aVec,aVec*mLin,'k--','linewidth',2);

hold on
co=get(gca,'colororder');
for kk=1:4
    ps(kk)=plot(aVec,funcs{kk}(aVec),'-','linewidth',2,...
        'color',co(kk,:));
    hold on
end

% legStr={'bound','first unbound','linear'};
% legend([ps(1:2) pFree],legStr,'location','southwest');

ylim([0 .5]);
xlim([0 .6]);

set(gca,'xgrid','on','ygrid','on','fontsize',10,'box','on',...
    'linewidth',1);

xlabel('scattering length $(a_\mathrm{HO})$','interpreter','latex');
ylabel('binding energy ($\hbar \omega$)','interpreter','latex');

%%

subplot(223)
p7=plot(Bvec,feshbach_97(Bvec),'-','LineWidth',2); 
hold on
p5=plot(Bvec,feshbach_95(Bvec),'-','LineWidth',2); 

ylim([-200 400]);

xlim([203 222]);

s7 = ['$|ab\rangle : 166.978a_0,~6.910~\mathrm{G},~202.15~\mathrm{G}$'];
s5 = ['$|ac\rangle : 174a_0,~9.7~\mathrm{G},~224.21\mathrm{G}$'];

sR = ['$a(B) = a_\mathrm{bg}\left(1-\frac{\Delta}{B-B_0}\right)$'];


legend([p7,p5],{s7,s5},'interpreter','latex',...
    'fontsize',9,'location','southeast');

text(.02,.98,sR,'interpreter','latex','units','normalized','fontsize',12,...
    'verticalalignment','top');

set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');

sL='PhysRevLett.90.230404, arxiv 2101.02700';
t=uicontrol('style','text','string',sL,'backgroundcolor','w',...
    'horizontalalignment','left','fontsize',8);
t.Position=[2 hF.Position(4)-15 300 15];

xlabel('magnetic field (G)');
ylabel('scattering length (a_0)');

drawnow;


%%
subplot(224);
co=get(gca,'colororder');
plot(Bvec,U_97*1e-3,'-','linewidth',2,'color',co(1,:))
hold on
plot(Bvec,U_97_linear*1e-3,'--','linewidth',2,'color',co(1,:))


plot(Bvec,U_95*1e-3,'-','linewidth',2,'color',co(2,:))
plot(Bvec,U_95_linear*1e-3,'--','linewidth',2,'color',co(2,:))

xlim([203 222]);
ylim([-75 75]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontname','times','fontsize',12);

xlabel('magnetic field (G)');
ylabel('s-wave shift (kHz)');

legend({'$U_{97}$','$U_{97}$ linear','$U_{95}$','$U_{95}$ linear'},'interpreter','latex',...
    'fontsize',9,'location','southeast');
%}

%%

data=load('G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\2021.09.29 Composite S-wave\analyzed_data.mat');

data2=struct;
data2.f_single = 44.379059736430;
data2.f_double= -31.82;
data2.B_fit = 199.31;

hF2=figure(2);
clf

set(gcf,'color','w');
hF2.Position=[1000 50 800 250];
hF2.Name = 'data_differential_swave';

subplot(121);
co=get(gca,'colororder');
plot(Bvec,U_97*1e-3,'-','linewidth',2,'color',co(1,:))
hold on
plot(Bvec,U_97_linear*1e-3,'--','linewidth',2,'color',co(1,:))


plot(Bvec,U_95*1e-3,'-','linewidth',2,'color',co(2,:))
plot(Bvec,U_95_linear*1e-3,'--','linewidth',2,'color',co(2,:))

xlim([203 222]);
ylim([-75 75]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontname','times','fontsize',10);

xlabel('magnetic field (G)');
ylabel('s-wave shift (kHz)');

legend({'$U_{97}$','$U_{97}$ linear','$U_{95}$','$U_{95}$ linear'},'interpreter','latex',...
    'fontsize',9,'location','southeast');


subplot(122);
co=get(gca,'colororder');
plot(Bvec,(U_95-U_97)*1e-3,'k-','linewidth',2)
hold on
plot(Bvec,(U_95_linear-U_97_linear)*1e-3,'k--','linewidth',2)

plot(data.B_fit,data.f_double,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8)
plot(data2.B_fit,data2.f_double,'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5,...
    'linewidth',2,'markersize',8)

xlim([195 222]);
ylim([-100 100]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontname','times','fontsize',10);

xlabel('magnetic field (G)');
ylabel('differential s-wave shift (kHz)');

legend({'$U_{97}-U_{95}$','$U_{97}-U_{95}$ linear'},'interpreter','latex',...
    'fontsize',9,'location','southeast');
end

