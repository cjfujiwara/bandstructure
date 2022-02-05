function differential_swave_shift4
out     = constants;
a0      = out.a0;
fr      = out.fr;
U       = 200;

% Band gap
% Delta   = sqrt(4*U)*fr
Delta   = [122529.554801461];

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

% Get the feshbach resonance
[a97,params_97] = feshbach_97(Bvec-dB);
[a95,params_95] = feshbach_95(Bvec-dB);
 
U_97 = Delta*int_func(a97/aho);
U_95 = Delta*int_func(a95/aho);

mLin = fout1.A;
mLin = sqrt(2/pi);

U_97_linear = Delta*(a97/aho)*mLin;
U_95_linear = Delta*(a95/aho)*mLin;

str = ['$\Delta = ' num2str(round(Delta*1e-3,2)) '~\mathrm{kHz}$' newline ...
    '$|ac\rangle : ' num2str(params_95.a_bg) 'a_0,~' num2str(params_95.B0) '~\mathrm{G},~' num2str(params_95.Delta) '~\mathrm{G}$' newline ... 
    '$|ab\rangle : ' num2str(params_97.a_bg) 'a_0,~' num2str(params_97.B0) '~\mathrm{G},~' num2str(params_97.Delta) '~\mathrm{G}$'];



%%

data=load('G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Composite S-wave RF\swave_RF_ALL.mat');
data=data.data_process;

data_raman=load('G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Composite S-wave RF\data_swave_raman.mat');
data_raman=data_raman.data_process;

hF2=figure;
clf

set(gcf,'color','w');
hF2.Position=[1000 50 600 600];
hF2.Name = 'data_differential_swave_RF_only';

subplot(211);
co=get(gca,'colororder');
plot(Bvec,(U_95-U_97)*1e-3,'k-','linewidth',2)
hold on
plot(Bvec,(U_95_linear-U_97_linear)*1e-3,'k--','linewidth',2)


errorbar(data.B,data.f1,...
    data.s1,data.s1,...
    data.B_err,data.B_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');



xlim([197 223]);
ylim([-100 100]);
ylim([-150 170]);

set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontsize',14);

xlabel('magnetic field (G)');
ylabel('differential s-wave shift (kHz)');

legend({'$U_{95}-U_{97}$','$U_{95}-U_{97}$ linear'},'interpreter','latex',...
    'fontsize',12,'location','southeast');

% legend({'$U_{95}-U_{97}$'},'interpreter','latex',...
%     'fontsize',12,'location','southeast');

subplot(212);
co=get(gca,'colororder');

% 97
plot(Bvec,U_97*1e-3,'-','linewidth',2,'color',co(1,:))
hold on
plot(Bvec,U_97_linear*1e-3,'--','linewidth',2,'color',co(1,:))


xlim([196 223]);
ylim([-135 135]);

set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontsize',12,...
    'fontname','helvetica');

xlabel('magnetic field (G)');
ylabel('s-wave shift (kHz)');


a95_data_ideal = feshbach_95(data.B-dB);
U_95_data_ideal = Delta*int_func(a95_data_ideal/aho);


% p1=plot(data.B,U_95_data_ideal*1e-3-data.f1,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
%     'linewidth',2,'markersize',8);

p1=errorbar(data.B,U_95_data_ideal*1e-3-data.f1,...
    data.s1,data.s1,...
    data.B_err,data.B_err,...
    'o','markerfacecolor',co(1,:),...
    'color',co(1,:)*.5,'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8); 


% p2=plot(data_raman.B,data_raman.f1,'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5,...
%     'linewidth',2,'markersize',8);


p2=errorbar(data_raman.B,data_raman.f1,...
    data_raman.s1,data_raman.s1,...
    data_raman.B_err,data_raman.B_err,...
    'o','markerfacecolor',co(2,:),...
    'color',co(2,:)*.5,'markeredgecolor',co(2,:)*.5,...
    'linewidth',2,'markersize',8); 

legend([p1 p2],{'RF','Raman'},'interpreter','latex',...
    'fontsize',12,'location','southeast');

text(.98,.98,str,'units','normalized','interpreter','latex','verticalalignment','top','horizontalalignment','right');


end

