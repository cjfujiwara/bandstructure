function doit

%% k1,k2 Energy Map
eng=@(k1,k2) -2*cos(pi*k1)-2*cos(pi*k2);    % Energy Functional
k1V = linspace(-3,3,1e3+1);                 % k1 vector
k2V = k1V;                                  % k2 vector
[kk1,kk2]=meshgrid(k1V,k2V);                % meshgrid of k1 and k2
engMap = eng(kk1,kk2);                      % image of energy
klim =[-2 2];                               % some plotting limits
%% Brillouin Zones
bz1_k1k2 = [1 1 -1 -1; -1 1 1 -1];          % First Brillouin Zone
bz2_k1k2 = [2 0 -2 0; 0 2 0 -2];            % Second Brillouin Zone
%% q,P Energy Map
% q:= k1-k2, P:=k1+k2

eng_qP=@(q,P) -4*cos(pi*q/2).*cos(pi*P/2);  % Energy Functional
qV = linspace(-6,6,1e3+1);                  % q vector (relative)
PV = qV;                                    % P vector (total)
[qq,PP]=meshgrid(qV,PV);                    % meshgrid of p and Q
engMap_qP = eng_qP(qq,PP);                  % image of energy

%% Change of Basis Vectors
% Unitary matrix to transform k1,k2, to q,P and visa-versa
U = [1 -1;1 1];
%% Choose Vectors

% Initial State |k1,k2>
k1k2 = [-.8;-.5];

% qP State |q,P>
qP_1=U*k1k2;

% qP+- State |q,P+-2> from umklapp scattering
qP_Plus = qP_1 + [0;2];
qP_Nega = qP_1 - [0;2];

% k3k4 States from Umklapp
k3k4_Plus = inv(U)*qP_Plus;
k3k4_Nega = inv(U)*qP_Nega;

%% Make Default Pictures
hF=figure(20);
clf
hF.Color='w';

baseMap = gray(256);        % Or any built-in colormap
white = ones(1, 3);           % RGB for white
alpha = 0.5;                  % 0 = white, 1 = full color
mutedMap = alpha * baseMap + (1 - alpha) * white;
colormap(mutedMap);


ax1 = subplot(121);
hImg = imagesc(k1V,k2V,engMap);
hold on
levels=linspace(-4,4,11);
contour(k1V,k2V,engMap,levels,'-','color',[.2 .2 .2])
cc=colorbar;
xlabel('$k_1~ [\pi/a]$','interpreter','latex');
ylabel('$k_2 ~[\pi/a]$','interpreter','latex');
cc.Label.String='energy [t] : $\mathcal{E} = -2 t \cos(k_1) -2t\cos(k_2)$';
cc.Label.Interpreter='latex';
plot(polyshape([1 0]*bz1_k1k2,[0 1]*bz1_k1k2),'facealpha',0,'linewidth',2);
plot(polyshape([1 0]*bz2_k1k2,[0 1]*bz2_k1k2),'facealpha',0,'linewidth',2,'linestyle','--');
axis equal tight
xlim(klim)
ylim(klim);
set(gca,'ydir','normal');
set(gca,'fontsize',16)
set(ax1,'fontsize',16)
set(ax1,'YTick',[-2 -1 0 1 2],'YTick',-2:1:2)
set(ax1,'YGrid','on','xgrid','on','fontname','times')
 co=get(gca,'colororder');
 title('$\left|k_1,k_2\right\rangle$-space','interpreter','latex')

ax2 = subplot(122);
hImg = imagesc(qV,PV,engMap_qP);
hold on
levels=linspace(-4,4,11);
contour(qV,PV,engMap_qP,levels,'-','color',[.2 .2 .2])
cc=colorbar;
xlabel('$q:=k_1-k_2 ~[\pi/a]$','interpreter','latex');
ylabel('$P:=k_1+k_2~ [\pi/a]$','interpreter','latex');
cc.Label.String='energy [t] : $\mathcal{E} = -4 t \cos(q/2)\cos(P/2)$';

cc.Label.Interpreter='latex';
axis equal tight
xlim([-2.5 2.5])
ylim([-2.5 2.5]);
set(gca,'ydir','normal');
plot(polyshape([1 0]*U*bz1_k1k2,[0 1]*U*bz1_k1k2),'facealpha',0,'linewidth',2,'linestyle','-');
plot(polyshape([1 0]*U*bz2_k1k2,[0 1]*U*bz2_k1k2),'facealpha',0,'linewidth',2,'linestyle','--');
 title('$\left|q,P\right\rangle$-space','interpreter','latex')

set(gca,'fontsize',16)
set(gca,'YTick',[-2 -1 0 1 2],'YTick',-2:1:2)
set(gca,'YGrid','on','xgrid','on','fontname','times')
 %% Specific Vectors



 % Draw states that converse total momentum and energy (but not relative).
k1_dummy = linspace(-1,1,100);
P0 =1;
k2_dummy = P0-k1_dummy;
k1mod = mod(k1_dummy+1,2)-1;
k2mod = mod(k2_dummy+1,2)-1;
ds=diff(k2mod);
[~,ind]=max(abs(ds));
plot(k1mod(1:ind),k2mod(1:ind),'-','parent',ax1,'linewidth',2,...
    'color',co(3,:));
plot(k1mod(ind+1:end-1),k2mod(ind+1:end-1),'-','parent',ax1,'linewidth',2,....
    'color',co(3,:));

 % Draw states that converse total momentum (but not relative).
% This includes all possible umklapp events
k1_dummy = linspace(-1,1,100);
P0 = k1k2(1)+k1k2(2);
k2_dummy = P0-k1_dummy;
k1mod = mod(k1_dummy+1,2)-1;
k2mod = mod(k2_dummy+1,2)-1;
ds=diff(k2mod);
[~,ind]=max(abs(ds));
plot(k1mod(1:ind),k2mod(1:ind),'-','parent',ax1,'linewidth',2,...
    'color',.8*co(1,:));
plot(k1mod(ind+1:end-1),k2mod(ind+1:end-1),'-','parent',ax1,'linewidth',2,....
    'color',.8*co(2,:));


% Initial Momentum (k1,k2);
plot(k1k2(1),k1k2(2),'ko','markerfacecolor',co(1,:),...
    'linewidth',1,'color','k','parent',ax1,...
    'markersize',6,'markersize',10);


% Umklapp Processes in k1,k2
plot(k3k4_Plus(1),k3k4_Plus(2),'marker','s','color','k','parent',ax1,...
    'markerfacecolor',co(2,:),'markersize',8)
plot(k3k4_Nega(1),k3k4_Nega(2),'marker','s','color','k','parent',ax1,...
    'markerfacecolor',co(2,:),'markersize',8)
%%

% Special Umklapp that conserve energy
plot([-1 1],[1 1],'-','parent',ax2,'linewidth',2,...
    'color',co(3,:));
plot([-1 1],[-1 -1],'-','parent',ax2,'linewidth',2,...
    'color',co(3,:));


qp_Dummy1=U*[k1mod(1) k1mod(ind);k2mod(1) k2mod(ind)];
qp_Dummy2=U*[k1mod(ind+1) k1mod(end-1);k2mod(ind+1) k2mod(end-1)];

 plot([1 0]*qp_Dummy1,[0 1]*qp_Dummy1,'-','parent',ax2,'linewidth',2,...
    'color',.8*co(1,:));
 plot([1 0]*qp_Dummy2,[0 1]*qp_Dummy2,'-','parent',ax2,'linewidth',2,...
    'color',.8*co(2,:));


% Umklapp Plus 1
plot([-1 1],[1 1],'-','parent',ax2,'linewidth',2,...
    'color',co(3,:));
plot([-1 1],[-1 -1],'-','parent',ax2,'linewidth',2,...
    'color',co(3,:));

% Initial Momentum (k1,k2);
plot(qP_1(1),qP_1(2),'ko','markerfacecolor',co(1,:),...
    'linewidth',1,'color','k','parent',ax2,...
    'markersize',10);

% Umklapp Processes in qP
plot(qP_Plus(1),qP_Plus(2),'marker','s','color','k','parent',ax2,...
    'markerfacecolor',co(2,:),'markersize',8)
plot(qP_Nega(1),qP_Nega(2),'marker','s','color','k','parent',ax2,...
    'markerfacecolor',co(2,:),'markersize',8)


end

