function differential_swave_shift2

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
a97 = feshbach_97(Bvec-dB);
a95 = feshbach_95(Bvec-dB);
 
U_97 = Delta*int_func(a97/aho);
U_95 = Delta*int_func(a95/aho);

mLin = fout1.A;
mLin = sqrt(2/pi);

U_97_linear = Delta*(a97/aho)*mLin;
U_95_linear = Delta*(a95/aho)*mLin;


%%

data=load('G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\2021.09.29 Composite S-wave\analyzed_data.mat');

data2=struct;
data2.f_single = 46612531;
data2.f_double= -31.82;
data2.B_fit = 198.56;

% 08/12 R10, 08/13 R04, 08/19 07,08, 08/20 03
data3=struct;
data3.f_single = [47041535 47128978 46871404 47162160 46958957];
data3.f_double= [-71.68 -94.27 -49 -98.11 -58.53];
data3.B_fit = [201.07 201.59 200.08 201.784 200.591];

% Raman
data4 = struct;
f1 = [-133.7 -133.8 -133.2 -132.5 -133 -133.3 -133.1];
f2 = [-73.7 -45.8 -68.7 -152.7 -169.4 -163.2 -154.2];
B = [198 200 198.5 205.5 204 204.5 205]+.11;
data4.f_double = f2-f1;
data4.B_fit = B;

hF2=figure(2);
clf

set(gcf,'color','w');
hF2.Position=[1000 50 800 600];
hF2.Name = 'data_differential_swave';

subplot(211);
co=get(gca,'colororder');
plot(Bvec,(U_95-U_97)*1e-3,'k-','linewidth',2)
hold on
plot(Bvec,(U_95_linear-U_97_linear)*1e-3,'k--','linewidth',2)

plot(data.B_fit,data.f_double,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8)
plot(data2.B_fit,data2.f_double,'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5,...
    'linewidth',2,'markersize',8)
plot(data3.B_fit,data3.f_double,'o','markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.5,...
    'linewidth',2,'markersize',8)

xlim([195 222]);
ylim([-100 100]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontname','times','fontsize',10);

xlabel('magnetic field (G)');
ylabel('differential s-wave shift (kHz)');

legend({'$U_{97}-U_{95}$','$U_{97}-U_{95}$ linear'},'interpreter','latex',...
    'fontsize',9,'location','southeast');


subplot(212);
co=get(gca,'colororder');

% 97
plot(Bvec,U_97*1e-3,'-','linewidth',2,'color',co(1,:))
hold on
plot(Bvec,U_97_linear*1e-3,'--','linewidth',2,'color',co(1,:))

% 95
plot(Bvec,U_95*1e-3,'-','linewidth',2,'color',co(2,:))
plot(Bvec,U_95_linear*1e-3,'--','linewidth',2,'color',co(2,:))

xlim([195 222]);
ylim([-75 140]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontname','times','fontsize',10);

xlabel('magnetic field (G)');
ylabel('s-wave shift (kHz)');


a95_data_ideal = feshbach_95(data.B_fit-dB);
U_95_data_ideal = Delta*int_func(a95_data_ideal/aho);

a95_data_ideal2 = feshbach_95(data2.B_fit-dB);
U_95_data_ideal2 = Delta*int_func(a95_data_ideal2/aho);

a95_data_ideal3 = feshbach_95(data3.B_fit-dB);
U_95_data_ideal3 = Delta*int_func(a95_data_ideal3/aho);

plot(data.B_fit,U_95_data_ideal*1e-3-data.f_double,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8)

plot(data2.B_fit,U_95_data_ideal2*1e-3-data2.f_double,'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5,...
    'linewidth',2,'markersize',8)

plot(data3.B_fit,U_95_data_ideal3*1e-3-data3.f_double,'o','markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.5,...
    'linewidth',2,'markersize',8)

plot(data4.B_fit,data4.f_double,'o','markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5,...
    'linewidth',2,'markersize',8)

legend({'$U_{97}$','$U_{97}$ linear','$U_{95}$','$U_{95}$ linear'},'interpreter','latex',...
    'fontsize',9,'location','southeast');



end

