function showTunnelingDepth(npt)

fme = figure;
clf
fme.Color='w';


X = [npt.depth];

t_11 = npt.Tunneling(1,1,:);t_11=t_11(:);
t_12 = npt.Tunneling(1,2,:);t_12=t_12(:);
t_13 = npt.Tunneling(1,3,:);t_13=t_13(:);

t_21 = npt.Tunneling(2,1,:);t_21=t_21(:);
t_22 = npt.Tunneling(2,2,:);t_22=t_22(:);
t_23 = npt.Tunneling(2,3,:);t_23=t_23(:);

t_31 = npt.Tunneling(3,1,:);t_31=t_31(:);
t_32 = npt.Tunneling(3,2,:);t_32=t_32(:);
t_33 = npt.Tunneling(3,3,:);t_33=t_33(:);

co = get(gca,'colororder');

subplot(311);
plot(X,npt.fr*t_11,'-','linewidth',1,'color',co(1,:));
hold on
plot(X,npt.fr*t_12,'--','linewidth',1,'color',co(1,:));
plot(X,npt.fr*t_13,'-.','linewidth',1,'color',co(1,:));
ylabel('tunneling rate (Hz)','interpreter','latex');
xlabel('depth ($E_R$)','interpreter','latex')
% xlim([1 10]);
title('$s$-band','interpreter','latex')

subplot(312);
plot(X,npt.fr*t_21,'-','linewidth',1,'color',co(2,:));
hold on
plot(X,npt.fr*t_22,'--','linewidth',1,'color',co(2,:));
plot(X,npt.fr*t_23,'-.','linewidth',1,'color',co(2,:));
ylabel('tunneling rate (Hz)','interpreter','latex');
xlabel('depth ($E_R$)','interpreter','latex')
% xlim([1 10]);
title('$p$-band','interpreter','latex')

subplot(313);
plot(X,npt.fr*t_31,'-','linewidth',1,'color',co(3,:));
hold on
plot(X,npt.fr*t_32,'--','linewidth',1,'color',co(3,:));
plot(X,npt.fr*t_33,'-.','linewidth',1,'color',co(3,:));
ylabel('tunneling rate (Hz)','interpreter','latex');
xlabel('depth ($E_R$)','interpreter','latex')
% xlim([1 10]);
title('$d$-band','interpreter','latex')
end

