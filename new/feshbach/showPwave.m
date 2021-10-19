function hF = showPwave
%% Resonance

Bvec=linspace(195,202,1000);
hF=figure(231);
set(hF,'color','w','name','feshbach_resonance_77');
hF.Position=[575 455 575 250];
clf

[V0,V1]=feshbach_77(Bvec);

p7=plot(Bvec,V0*1e-6,'-','LineWidth',2); 
hold on
p5=plot(Bvec,V1*1e-6,'-','LineWidth',2); 

ylim([-2e2 2e2]);

xlim([197 200]);

s7 = ['$M_L=0 : -(108.0 a_0)^3,~-19.89~\mathrm{G},~198.803\mathrm{G}$'];
s5 = ['$M_L=1 : -(107.35 a_0)^3,~-19.54~\mathrm{G},~198.300\mathrm{G}$'];
sR = ['$V(B) = V_\mathrm{bg}\left(1-\frac{\Delta}{B-B_0}\right)$'];


legend([p7,p5],{s7,s5},'interpreter','latex',...
    'fontsize',9,'location','southeast');

text(.02,.98,sR,'interpreter','latex','units','normalized','fontsize',12,...
    'verticalalignment','top');

set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');

sL='arxiv 2101.02700';
t=uicontrol('style','text','string',sL,'backgroundcolor','w',...
    'horizontalalignment','left','fontsize',8);
t.Position=[2 2 100 15];

xlabel('magnetic field (G)');
ylabel('scattering length (100 a_0)^3');

drawnow;

end

