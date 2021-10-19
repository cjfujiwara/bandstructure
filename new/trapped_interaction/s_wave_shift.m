function s_wave_shift

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

% Get the feshbach resonance
a97 = feshbach_97(Bvec+dB);
a95 = feshbach_95(Bvec+dB);
 
U_97 = Delta*int_func(a97/aho);
U_95 = Delta*int_func(a95/aho);

mLin = fout1.A;
mLin = sqrt(2/pi);

U_97_linear = Delta*(a97/aho)*mLin;
U_95_linear = Delta*(a95/aho)*mLin;


%%

data=load('G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\2021.09.29 Composite S-wave\analyzed_data.mat');

hF2=figure(3);
clf

set(gcf,'color','w');
hF2.Position=[1150 300 750 250];

co=get(gca,'colororder');
hold on

a95_data_ideal = feshbach_95(data.B_fit+dB);
U_95_data_ideal = Delta*int_func(a95_data_ideal/aho);

plot(Bvec,U_97*1e-3,'-','linewidth',2,'color',co(1,:))
plot(Bvec,U_97_linear*1e-3,'--','linewidth',2,'color',co(1,:))

plot(data.B_fit,U_95_data_ideal*1e-3-data.f_double,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8)

xlim([203 222]);ylim([-80 20]);
set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on','fontname','times','fontsize',10);

xlabel('magnetic field (G)');
ylabel('97 s-wave shift (kHz)');

legend({'sat.','linear','data'},'interpreter','latex',...
    'fontsize',9,'location','southeast');
end

