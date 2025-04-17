function doit

%% k1,k2 Energy Map

eng=@(k1,k2) -2*cos(pi*k1)-2*cos(pi*k2);
k1V = linspace(-3,3,1e3+1);
k2V = k1V;
[kk1,kk2]=meshgrid(k1V,k2V);
engMap = eng(kk1,kk2);
klim =[-3 3];


%% q,P Energy Map

eng_qP=@(q,P) -4*cos(pi*q/2).*cos(pi*P/2);

qV = linspace(-6,6,1e3+1);
PV = qV;
[qq,PP]=meshgrid(qV,PV);
engMap_qP = eng_qP(qq,PP);


% Convert (k1,k2) to (q,P) (q=k1-k2,P=k1+k2);

U_kk_2_qP = [1 -1;1 1];
U_qP_2_kk = inv(U_kk_2_qP);

k1k2_2_qP =@(k1,k2) [1 -1;1 1]*[k1;k2];

klim_2_qp_a = k1k2_2_qP(klim(2),klim(2));
klim_2_qp_b = k1k2_2_qP(klim(1),klim(2));
klim_2_qp_c = k1k2_2_qP(klim(1),klim(1));
klim_2_qp_d = k1k2_2_qP(klim(2),klim(1));


klim_in_QP_x = [klim_2_qp_a(1) klim_2_qp_b(1) klim_2_qp_c(1) klim_2_qp_d(1)];
klim_in_QP_y = [klim_2_qp_a(2) klim_2_qp_b(2) klim_2_qp_c(2) klim_2_qp_d(2)];

kfbz_2_qp_a = k1k2_2_qP(1,1);
kfbz_2_qp_b = k1k2_2_qP(-1,1);
kfbz_2_qp_c = k1k2_2_qP(-1,-1);
kfbz_2_qp_d = k1k2_2_qP(1,-1);
kfbz_in_QP_x = [kfbz_2_qp_a(1) kfbz_2_qp_b(1) kfbz_2_qp_c(1) kfbz_2_qp_d(1)];
kfbz_in_QP_y = [kfbz_2_qp_a(2) kfbz_2_qp_b(2) kfbz_2_qp_c(2) kfbz_2_qp_d(2)];

%% Choose Vectors
% Initial Momentum
k1 = .7;
k2 = -0.5;

str1=['$(' num2str(k1) ',' num2str(k2) ')$'];

% qP Initial
v1=U_kk_2_qP*[k1;k2];
q1 = v1(1);P1=v1(2);
disp(v1);

% Final possible umklapp
q2 = q1;
P2_p = P1+2;
P2_n = P1-2;

disp('k1k2')
disp(eng(k1,k2))

disp('q1P1')
disp(eng_qP(q1,P1))

disp(eng_qP(q1,P1+2))

% Final possible umklapp in k1,k2
% % qP Final
v2=U_qP_2_kk*[q2;P2_n];
v3=U_qP_2_kk*[q2;P2_p];

%% Check Energies

% Find q prime that corresponds to an umklapp event and conserves total
% final energy
E1 = eng_qP(q1,P1);
n=4;
qprime = fzero(@(qprime) eng_qP(qprime,P1-n)-E1,0);

k_umklapp = U_qP_2_kk*[qprime;P1-n];

while abs(k_umklapp(1))>1
    k_umklapp(1)=k_umklapp(1)-2*sign(k_umklapp(1));
end

while abs(k_umklapp(2))>1
    k_umklapp(2)=k_umklapp(2)-2*sign(k_umklapp(2));
end


% ind=find(E1-eng_qP(qprime,P1-2)<0,1)
% keyboard
% disp(eng(k1,k2))


% disp(eng(v2(1),v2(2)))

%% Make Default Pictures
hF=figure(20);
clf
hF.Color='w';
colormap(parula)


ax1 = subplot(121);
hImg = imagesc(k1V,k2V,engMap,'alphadata',.5);
hold on
levels=linspace(-4,4,10);
contour(k1V,k2V,engMap,levels,'k-')
cc=colorbar;
xlabel('k_1');
ylabel('k_2');
cc.Label.String='energy (t)';
bz_1 = plot(polyshape([1 1 -1 -1],[-1 1  1 -1]),'facealpha',0,'linewidth',2);
bz_2 = plot(polyshape([2 0 -2 0],[0 2 0 -2]),'facealpha',0,'linewidth',2);
axis equal tight
xlim(klim)
ylim(klim);
set(gca,'ydir','normal');
set(gca,'fontsize',16)
set(ax1,'fontsize',16)
set(ax1,'YTick',[-2 -1 0 1 2],'YTick',-2:1:2)
set(ax1,'YGrid','on','xgrid','on')

% k=45;
% set(ax1, 'CameraUpVector', [sind(k), cosd(k), 0]); %% Specific Vectors


ax2 = subplot(122);
hImg = imagesc(qV,PV,engMap_qP,'alphadata',.5);
hold on
levels=linspace(-4,4,10);
contour(qV,PV,engMap_qP,levels,'k-')
cc=colorbar;
xlabel('$q=k_1-k_2$','interpreter','latex');
ylabel('$P=k_1+k_2$','interpreter','latex');
cc.Label.String='energy (t)';

axis equal tight
xlim([-4 4])
ylim([-4 4]);
set(gca,'ydir','normal');

k_limits = plot(polyshape(klim_in_QP_x,klim_in_QP_y),'facealpha',0,'linewidth',1,'linestyle','--');
k_fbz = plot(polyshape(kfbz_in_QP_x,kfbz_in_QP_y),'facealpha',0,'linewidth',1,'linestyle','-');

set(gca,'fontsize',16)
set(gca,'YTick',[-2 -1 0 1 2],'YTick',-2:1:2)
set(gca,'YGrid','on','xgrid','on')
 %% Specific Vectors

% Initial Momentum (k1,k2);
plot([0 k1],[0 k2],'linewidth',2,'color','k','parent',ax1);

% Initial Momentum (q,P)
plot([0 q1],[0 P1],'-','linewidth',2','color','k','parent',ax2)

% Umklapp Processes
plot([0 q2],[0 P2_n],'-','linewidth',2,'color','r','parent',ax2)
plot([0 q2],[0 P2_p],'-','linewidth',2,'color','g','parent',ax2)


% Umklapp Processes in k1,k2
plot([0 v2(1)],[0 v2(2)],'-','linewidth',2,'color','r','parent',ax1)
plot([0 v3(1)],[0 v3(2)],'-','linewidth',2,'color','g','parent',ax1)
end

