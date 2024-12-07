
Omega=2;
t = 563;


N=601;
nc = (N+1)/2;

% Change of basis matrix
u1 = eye(N);
u2 = flip(u1,1);
u1(nc:end,nc:end)=-u1(nc:end,nc:end);
U=u1+u2;
U=U/sqrt(2);
U(nc,nc)=1;

% Potential energy matrix
x=[-(nc-1):1:(nc-1)];
X2=diag(x.^2);

% Curvature Matrix
T=makeTmatrix(N,1);

Tstar = U*T*ctranspose(U);

Tstar(abs(Tstar)<1e-3)=0;

%% Construct hamiltonian

% Hamiltonian Normal
H1 = T*t+Omega*X2;

% Hamiltonian in new basis
H2 = U*H1*ctranspose(U);
H2(abs(H2)<1e-3)=0;

H2a = H2(1:nc,1:nc);
H2b = H2((nc+1):end,(nc+1):end);
%% Compare Methods

% Original
[c1,d1]=eig(H1);
e1 = diag(d1);

% New Basis
[c2,d2]=eig(H2);
e2 = diag(d2);

% New Basis Even
[c2a,d2a]=eig(H2a);
e2a = diag(d2a);

% New Basis Odd
[c2b,d2b]=eig(H2b);
e2b = diag(d2b);

% Convert into original basis
C = U*blkdiag(c2a,c2b);

[e2c,inds] = sort([e2a; e2b],'ascend');

C = C(:,inds);

hF=figure(20);
hF.Color='w';
clf
co=get(gca,'colororder');
plot(e1,'ko-','markerfacecolor',co(1,:),'markersize',6);
hold on
plot(e2,'.-','markerfacecolor',co(2,:),'markersize',6);
plot(e2c,'ko','markerfacecolor',co(3,:),'markersize',4);

xlim([0 50])

figure(1);
clf
imagesc(x,1:N,abs(C)'.^2); set(gca,'YDir','normal')
xlabel('site');
ylabel('eigenindex');
cc=colorbar;
cc.Label.String='abs(\psi(n))^2';
title('2.5 Er + 60 Hz spectrum')
caxis([0 .05]);
xlim([-60 60])
ylim([1 100]);

%% Wavefunction more

w=real(npt.Wannier_X(:,1));
w=w/norm(w);

xw=npt.X_extended;

dx=xw(2)-xw(1);
dn = round(1/dx);
L = 10*dn;

xc = (length(xw)+1)/2;

i1 = xc-L;
i2 = xc+L;

w_sub = w(i1:i2);
w_sub = w_sub/norm(w_sub);
x_sub  = xw(i1:i2);
ind=10;
x_fine = linspace(x(1),x(end),round((x(end)-x(1))/dx)+1);
y_fine = zeros(length(x_fine),1);

psi_all = zeros(length(x_fine),size(C,2));
y_fine = zeros(length(x_fine),1);

tic
for nn=1:size(C,2)
    y_fine=y_fine*0;
    y_fine(1:dn:end)=C(:,nn);
    psi_all(:,nn) = conv(y_fine,w_sub,'same');
end
toc
%%
figure(3);
clf
imagesc(x_fine,1:N,abs(psi_all)'.^2); set(gca,'YDir','normal');

xlabel('site');
ylabel('eigenindex');
cc=colorbar;
cc.Label.String='abs(\psi(n))^2';
title('2.5 Er + 60 Hz spectrum')
caxis([0 1e-4]);
xlim([-60 60])
ylim([1 100]);
