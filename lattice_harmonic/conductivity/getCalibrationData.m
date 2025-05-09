function [outputArg1,outputArg2] = getCalibrationData

%% Load Data


trap_freq{1}=load('X:\Data\2024\2024.11\11.25\18 lattice xdt trap freq 2.5 ER req, 195 G (195,150) mW, xdtB spin pol, 2 V modulation 120 kHz ps amp try 3\figures\ixon_gaussdata.mat');
trap_freq{2}=load('X:\Data\2024\2024.12\12.17\04 ixon_Trap freq x-dir, ODT Powers 195,150 mW2.5,2.5,2.5 Er201.1G,65 mW 2 amp, 50 ms ramp\figures\ixon_gaussdata.mat');
trap_freq{3}=load('X:\Data\2024\2024.12\12.18\02 ixon_Trap freq, ODTs (195,150) mW, (2.5,2.5,2.5) Er, 201.1G, 65 mW, 2 V amp, 50 ms ramp\figures\ixon_gaussdata.mat');
trap_freq{4}=load('X:\Data\2025\2025.01\01.22\07 ixon_Trap freq, ODTs (195,150) mW, (2.5,2.5,2.5) Er, 201.1 G, 64.5 mW, 3 V amp, 50 ms ramp\figures\ixon_gaussdata.mat');
trap_freq{5}=load('X:\Data\2025\2025.02\02.11\13 ixon_Trap freq, ODTs (195,150) mW, (2.5,2.5,2.5) Er, 200 G, 65 mW, 3 V amp, 50 ms ramp, Vert Disp 5 V\figures\ixon_gaussdata.mat');
trap_freq{6}=load('X:\Data\2025\2025.03\03.06\11 Xdir ixon_Trap freq, ODTs (195,150) mW, (2.5,2.5,2.5) Er, 201.1 G, 54 mW, 3 V amp, 50 ms ramp\figures\ixon_gaussdata.mat');
trap_freq{7}=load('X:\Data\2025\2025.03\03.11\03 ixon_XDT and Lattice X Trap Freq, ODTs (198,88) mW, (2.5,2.5,2.5) Er, 201.1 G, 54 mW, 3 V amp, 50 ms ramp\figures\ixon_gaussdata.mat');

for kk=1:length(trap_freq)
   trap_freq{kk}=trap_freq{kk}.ixon_gaussdata; 
end

%%

decayFit = fittype('A*cos(2*pi*f*(t-t0))*exp(-(t-t0)/tau)+B',...
    'coefficients',{'A','t0','f','B','tau'},'independent','t');
fitopt = fitoptions(decayFit);


%%
hF=figure(20);
hF.Color='w';
clf
co=jet(length(trap_freq));
tt=linspace(0,40,1e3);

tVec = [];
fVec= zeros(length(trap_freq),2);
for kk=1:length(trap_freq)
    P=[trap_freq{kk}.Params];
    t=mean([P.ExecutionDate]);
    tstr=datestr(t,'YYYY-mm-dd');    
    tVec(kk)=t;
    
    subplot(length(trap_freq),2,2*(kk-1)+1);
    Y=[trap_freq{kk}.Xc]-trap_freq{kk}.Xc(1);Y=Y(:);
    X=[trap_freq{kk}.X];X=X(:);   
    
    
    plot(X,Y,'o','markerfacecolor',co(kk,:),...
        'markeredgecolor','k');
    title(tstr);

    hold on
    xlim([0 40]);
    ylabel('relative position (px)')
    ylim([-5 40]);
    
    Ag=-(max(Y)-min(Y))*0.5;
    t0g= -4;
    fg = 55*1e-3;
    Bg = mean(Y);
    taug = 40;
    fitopt.StartPoint=[Ag t0g fg Bg taug];
    fitopt.Weights = ones(length(Y),1);

    fitAll=fit(X,Y,decayFit,fitopt);
    pFAll=plot(tt,feval(fitAll,tt),'k--');    
    
    
    Xmin = 4;
    Xmax = 30;
    W=double([X>=Xmin].*[X<=Xmax]);W=W(:);
    
    fitopt.Weights=W;
    fitSome=fit(X,Y,decayFit,fitopt);
    pFSome=plot(tt,feval(fitSome,tt),'k-');

    cSome=confint(fitSome);
    dfSome=0.5*(cSome(2,3)-cSome(1,3));
    cAll=confint(fitAll);
    dfAll=0.5*(cAll(2,3)-cAll(1,3));
    
    legStr={[num2str(1e3*fitAll.f,'%.1f') ' \pm ' num2str(1e3*dfAll,'%.0f') ' Hz'],...
        [num2str(1e3*fitSome.f,'%.1f') ' \pm ' num2str(1e3*dfSome,'%.0f') ' Hz']};
    legend([pFAll pFSome],legStr,'orientation','horizontal','location','northeast');
    
    
    subplot(length(trap_freq),2,2*(kk-1)+2);
    Y2=[trap_freq{kk}.Xs];Y=Y(:);
    
    plot(X,Y2*16/80,'o','markerfacecolor',co(kk,:),...
        'markeredgecolor','k');
    ylabel('size (um)');
        xlim([0 40]);
        
    fVec(kk,1:2)=[fitSome.f dfSome];
    
    Xs(kk,1)=mean(Y2)*16/80;
    Xs(kk,2)=std(Y2)*16/80;


end
%     tVec=datetime(tVec,'convertfrom','datenum');
hF2=figure(21);
hF2.Color='w';
clf
cme=get(gca,'colororder');
errorbar(tVec,1e3*fVec(:,1),1e3*fVec(:,2),'ko','markerfacecolor',cme(1,:));
xlabel('calibration date');
ylabel('oscillation frequency (Hz)');
ylim([40 80]);
yyaxis right
errorbar(tVec,Xs(:,1),Xs(:,2),'ko','markerfacecolor',cme(2,:));
ylabel('second moment (um)');
ylim([5 10]);

datetick x
title('x lho oscillations');
% %% XDT
% load('X:\Data\2024\2024.12\12.03\ixon_ODT Powers 195,150 mW-0.5,-0.5,-0.5 Er201.1G,65 mW 2 amp, 50 ms ramp\figures\ixon_gaussdata.mat')
% load('X:\Data\2024\2024.12\12.17\05 XDT X-dir trap freq\figures\ixon_gaussdata.mat')
% load('X:\Data\2024\2024.12\12.18\03 ixon_Trap freq, ODTs (195,150) mW, (-0.5,-0.5,-0.5) Er, 201.1 G, 65 mW, 2 V amp, 50 ms ramp\figures\ixon_gaussdata.mat')
% load('X:\Data\2025\2025.01\01.13\03 ixon_Trap freq, ODTs (195,150) mW, (-0.5,-0.5,-0.5) Er, 201.1 G, 65 mW, 2 V amp, 50 ms ramp\figures\ixon_gaussdata.mat')

keyboard
end

