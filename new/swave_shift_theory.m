function swave_shift_theory
out     = constants;
a0      = out.a0;
fr      = out.fr;
U       = 200;

% Band gap
Delta_200   = sqrt(4*150)*fr;
Delta_100   = sqrt(4*100)*fr;
Delta_60   = sqrt(4*60)*fr;

% Harmonic oscillator length scale
aho=harmonic_length(U)/a0;

Bvec = linspace(190,225,1e4);

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

% Get the feshbach resonance
[a97,~] = feshbach_97(Bvec);
[a95,~] = feshbach_95(Bvec);
 
U_97 = int_func(a97/aho);
U_95 = int_func(a95/aho);



hF2=figure(11);
clf

set(gcf,'color','w');
hF2.Position=[500 100 800 400];
hF2.Name = 'data_differential_swave_RF_only';

subplot(121);
co=get(gca,'colororder');
plot(Bvec,Delta_200*(U_95-U_97)*1e-3,'-','linewidth',2,'color',co(1,:))
hold on
plot(Bvec,Delta_100*(U_95-U_97)*1e-3,'-','linewidth',2,'color',co(2,:))
plot(Bvec,Delta_60*(U_95-U_97)*1e-3,'-','linewidth',2,'color',co(3,:))
xlim([197 210]);
ylim([-100 100]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontsize',12,...
    'fontname','times');
xlabel('magnetic field (G)');
ylabel('$U_{95}-U_{97}$ (kHz) ','interpreter','latex');
legend({'200$E_R$','100$E_R$','60$E_R$'},'interpreter','latex',...
    'fontsize',12,'location','southeast');



subplot(122);
co=get(gca,'colororder');
plot(Bvec,Delta_200*(U_97)*1e-3,'-','linewidth',2,'color',co(1,:))
hold on
plot(Bvec,Delta_100*(U_97)*1e-3,'-','linewidth',2,'color',co(2,:))
plot(Bvec,Delta_60*(U_97)*1e-3,'-','linewidth',2,'color',co(3,:))
xlim([197 210]);
ylim([-100 100]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontsize',12,...
    'fontname','times');xlabel('magnetic field (G)');
legend({'200$E_R$','100$E_R$','60$E_R$'},'interpreter','latex',...
    'fontsize',12,'location','southeast');



end

