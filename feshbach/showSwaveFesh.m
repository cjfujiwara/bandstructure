function hF = showSwaveFesh
%% Resonance

Bvec=linspace(190,225,1000);
hF=figure(226);
set(hF,'color','w','name','feshbach_s_resonance');
hF.Position=[0 455 575 250];
clf
co=get(gca,'colororder');
p7=plot(Bvec,feshbach_97(Bvec),'-','LineWidth',2,'color',co(1,:)); 
hold on
p5=plot(Bvec,feshbach_95(Bvec),'-','LineWidth',2,'color',co(2,:)); 
p5b=plot(Bvec,feshbach_95_regal(Bvec),'--','LineWidth',2,'color',co(2,:)); 

ylim([-200 400]);

xlim([190 225]);

s7 = ['$|ab\rangle : 166.978a_0,~6.910~\mathrm{G},~202.15~\mathrm{G}$'];
s5 = ['$|ac\rangle : 167.3a_0,~7.2~\mathrm{G},~224.2\mathrm{G}$'];
s5b = ['$|ac\rangle : 174a_0,~9.7~\mathrm{G},~224.21\mathrm{G}$'];

sR = ['$a(B) = a_\mathrm{bg}\left(1-\frac{\Delta}{B-B_0}\right)$'];


legend([p7,p5,p5b],{s7,s5,s5b},'interpreter','latex',...
    'fontsize',9,'location','southeast');

text(.02,.98,sR,'interpreter','latex','units','normalized','fontsize',10,...
    'verticalalignment','top');

set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');

sL='PhysRevLett.90.230404, Ludewig, arxiv 2101.02700';
t=uicontrol('style','text','string',sL,'backgroundcolor','w',...
    'horizontalalignment','left','fontsize',8);
t.Position=[2 hF.Position(4)-15 300 15];

xlabel('magnetic field (G)');
ylabel('scattering length (a_0)');

drawnow;

end

