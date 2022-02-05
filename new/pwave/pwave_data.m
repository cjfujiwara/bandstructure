loc = 'G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Composite P-wave';

f100 =  'data_100Er.mat';
f200New = 'data_200Er_new.mat';
f200 =  'data_200Er.mat';
f300 =  'data_300Er.mat';

% f300 = 'data_300Er_old.mat';
d1 = load(fullfile(loc,f100));
d3 = load(fullfile(loc,f300));
d2New = load(fullfile(loc,f200New));
d2b= load(fullfile(loc,f200));



hf1=figure(10);
clf
hf1.Color='w';
hf1.Position = [100 100 800 400];
co=get(gca,'colororder');

p1=errorbar(d1.data_process.B,d1.data_process.df,...
    d1.data_process.df_err,d1.data_process.df_err,...
    d1.data_process.B_err,d1.data_process.B_err,...
    'o','markerfacecolor',co(1,:),...
    'markeredgecolor',.5*co(1,:),'color',co(1,:)*.5,...
    'linewidth',2,'markersize',8); 
hold on

p2=errorbar(d2.data_process.B,d2.data_process.df,...
    d2.data_process.df_err,d2.data_process.df_err,...
    d2.data_process.B_err,d2.data_process.B_err,...
    'o','markerfacecolor',co(2,:),...
    'markeredgecolor',.5*co(2,:),'color',co(2,:)*.5,...
    'linewidth',2,'markersize',8); 
hold on

p3=errorbar(d3.data_process.B,d3.data_process.df,...
    d3.data_process.df_err,d3.data_process.df_err,...
    d3.data_process.B_err,d3.data_process.B_err,...
    'o','markerfacecolor',co(3,:),...
    'markeredgecolor',.5*co(3,:),'color',co(3,:)*.5,...
    'linewidth',2,'markersize',8); 
hold on

xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');

ylim([-55 80]);
xlim([199 201.5]);
% xlim([round(min(all_B)-0.15,1) round(max(all_B)+0.15,1)]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',12,'xminorgrid','on','yminorgrid','on');

s1 = '$100~E_\mathrm{R}~(\Delta_\mathrm{O}=85.19~\mathrm{kHz})$';
s2 = '$200~E_\mathrm{R}~(\Delta_\mathrm{O}=122.54~\mathrm{kHz})$';
s3 = '$300~E_\mathrm{R}~(\Delta_\mathrm{O}=151.15~\mathrm{kHz})$';

legend([p1 p2 p3],{s1,s2,s3},'location','northwest',...
    'interpreter','latex','fontsize',12);
%%
hf2=figure(11);
hf2.Color='w';
hf2.Position(2:4) =[200 700 300];

clf
ind = 2;
src =d2New; 
co=get(gca,'colororder');

X = src.data_out(ind).X;
Y= src.data_out(ind).Y;
fout = src.data_out(ind).Fit;
mydir = src.data_out(ind).Directory
YStr = src.data_out(ind).YStr
hf2.Name=mydir;

% Find Unique Value    
[ux,ia,ib]=unique(X);    
Yu=zeros(length(ux),2);    
for kk=1:length(ux)
    inds=find(X==ux(kk));
    Yu(kk,1)=mean(Y(inds));
    Yu(kk,2)=std(Y(inds));       
end 

myco = co(2,:);
xx=linspace(-70,100,1e3);
plot(xx-fout.x1,feval(fout,xx),'-','linewidth',2,'color',myco*.5);
hold on
errorbar(ux-fout.x1,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',1,'markersize',8,'linewidth',2);    
hold on
set(gca,'xgrid','on','ygrid','on','fontname','times',...
    'fontsize',12);
xlim([-50 50]);
xlabel('$f-f_\mathrm{ZS}~(\mathrm{kHz})$','interpreter','latex','fontsize',16);
ylabel('$(N_\uparrow - N_\downarrow)/N$','interpreter','latex','fontsize',16);

%%
hf3=figure(12);
hf3.Color='w';
hf3.Position(2:4) =[200 700 300];
clf
ind = 2;
src =d3; 
co=get(gca,'colororder');

X = src.data_out(ind).X;
Y= src.data_out(ind).Y;
fout = src.data_out(ind).Fit;
mydir = src.data_out(ind).Directory
YStr = src.data_out(ind).YStr
hf3.Name=mydir;

% Find Unique Value    
[ux,ia,ib]=unique(X);    
Yu=zeros(length(ux),2);    
for kk=1:length(ux)
    inds=find(X==ux(kk));
    Yu(kk,1)=mean(Y(inds));
    Yu(kk,2)=std(Y(inds));       
end 

myco = co(3,:);
xx=linspace(-70,100,1e3);
plot(xx-fout.x1,feval(fout,xx),'-','linewidth',2,'color',myco*.5);
hold on
errorbar(ux-fout.x1,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
    'markeredgecolor',myco*.5,'color',myco,...
    'linewidth',1,'markersize',8,'linewidth',2);    
hold on
set(gca,'xgrid','on','ygrid','on','fontname','times',...
    'fontsize',12);
xlim([-50 50]);
xlabel('$f-f_\mathrm{ZS}~(\mathrm{kHz})$','interpreter','latex','fontsize',16);
ylabel('$(N_\uparrow - N_\downarrow)/N$','interpreter','latex','fontsize',16);