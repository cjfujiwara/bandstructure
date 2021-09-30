function [outputArg1,outputArg2] = differential_swave_shift(B,U)

out     = constants;
a0      = out.a0;
fr      = out.fr;
U       = 200;

% Band gap
Delta   = sqrt(4*U)*fr;

% Harmonic oscillator length scale
aho=harmonic_length(U)/a0;


Bvec = linspace(190,225,1e4);
aVec = linspace(-10,10,1e4);

% Calculat the resonance shift
[dB,~]=lattice_shift_79(U,U,U);

% Calculate the theoretical function
[funcs] = trapped_interaction;
int_func=funcs{2};

mLin = int_func(.05)/.05;

% Get the feshbach resonance
a97 = feshbach_97(Bvec+dB);
a95 = feshbach_95(Bvec+dB);
 
U_97 = Delta*int_func(a97/aho);
U_95 = Delta*int_func(a95/aho);

U_97_linear = Delta*(a97/aho)*sqrt(2/pi);
U_95_linear = Delta*(a95/aho)*sqrt(2/pi);

hF=figure(1);
clf

set(gcf,'color','w');
hF.Position=[50 50 800 800];

%% THeoery

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


%%

hF2=figure(2);
clf

set(gcf,'color','w');
hF2.Position=[600 50 800 350];

subplot(121);
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


subplot(122);
co=get(gca,'colororder');
plot(Bvec,(U_95-U_97)*1e-3,'k-','linewidth',2)
hold on
plot(Bvec,(U_95_linear-U_97_linear)*1e-3,'k--','linewidth',2)


xlim([203 222]);
ylim([25 100]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontname','times','fontsize',12);

xlabel('magnetic field (G)');
ylabel('differential s-wave shift (kHz)');

legend({'$U_{97}-U_{95}$','$U_{97}-U_{95}$ linear'},'interpreter','latex',...
    'fontsize',9,'location','north');
end
