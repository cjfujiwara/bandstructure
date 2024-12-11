x = pi*npt.X_extended';

wm = real(npt.Wannier_X(:,iBr,iU));        
wm = wm/sqrt(trapz(x,wm.*wm));

wn = real(npt.Wannier_X(:,iBc,iU));        
wn = wn/sqrt(trapz(x,wn.*wn));
% 
% 
% trapz(x,wm.*wm)
% trapz(x,wn.*wn)
% 
% trapz(x,wm.*x.*wn)

figure(29);
clf
plot(x,wm);
hold on
plot(x,wn);
xlim([-10 10]);



k = npt.K_extended;
dk = k(2)-k(1);

wmk = npt.Wannier_K(:,iBr,iU);
wmk = wmk/sqrt(trapz(k,conj(wmk).*wmk));
dwmk = gradient(wmk,dk);

wnk = npt.Wannier_K(:,iBc,iU);
wnk = wnk/sqrt(trapz(k,conj(wnk).*wnk));
dwnk = gradient(wnk,dk);       
                
%%
XL=[.99 1.01];

figure(30);
clf

subplot(131);
plot(k,real(wmk));
hold on;
plot(k,imag(wnk));
xlim([-2.5 2.5]);
xlim(XL);
xlabel('k');
ylabel('w(k)');

subplot(132);
plot(k,real(dwmk));
hold on;
plot(k,imag(dwnk));
xlim([-2.5 2.5]);
xlim(XL);
xlabel('k');
ylabel('w''(k)');

% Caclulate matrix element
nn=pi*1;
A = conj(real(wmk)).*imag(dwnk)*1i;
B = conj(real(dwmk)).*imag(wnk)*1i;
OSC = exp(-2*1i*k*nn);
C = (A-B)*0.5*1i.*OSC;
D=trapz(k,C);
disp(['dipole numeric ' num2str(real(D))]);


subplot(133);
plot(k,conj(wmk),'.-');
ylim([-.1 1]);
yyaxis right
plot(k,imag(dwnk),'.-');
xlim(XL);
ylim([- 100 2]);




%%



kt=linspace(-5,5,1e6+1)';
% kt=k;
dkt = kt(2)-kt(1);

% 
% wm_t = (sign(kt+1)-sign(kt-1))/(2*sqrt(2));
% wn_t = -1i*(-sign(kt-2)+sign(kt-1)-sign(kt+2)+sign(kt+1))/(2*sqrt(2));

foo=@(K) (heaviside(K+1)-heaviside(K-1))/(sqrt(2));



wm_t = (heaviside(kt+1)-heaviside(kt-1))/(sqrt(2));
wn_t = -1i*(-heaviside(kt-2)+heaviside(kt-1)-heaviside(kt+2)+heaviside(kt+1))/sqrt(2);

i1=find(kt==-1,1);
i2 = find(kt==1,1);

wm_t(i1) = 1/2;
wm_t(i2) = 1/2;

wn_t(i1) = 1i/2;
wn_t(i2) = -1i/2;


dwm_t = gradient(wm_t,dkt);
dwn_t = gradient(wn_t,dkt);

figure(31);
clf
subplot(131);
plot(kt,wm_t,'.-');
hold on
plot(kt,imag(wn_t));
xlim([-2.5 2.5]);
xlim(XL);

subplot(132);
plot(kt,dwm_t);
hold on;
plot(kt,imag(dwn_t));
xlim([-2.5 2.5]);
xlim(XL);


nn=pi*1;
At = conj(wm_t).*dwn_t;
Bt = conj(dwm_t).*wn_t;
OSCt = exp(-2*1i*kt*nn);
Ct = (At-Bt)*0.5*1i.*OSCt;
Dt = trapz(kt,Ct);
disp(Dt);

subplot(133);
plot(kt,conj(wm_t),'.-');
ylim([-.1 1]);
yyaxis right
plot(kt,imag(dwn_t),'.-');
ylim([- 100 2]);
xlim(XL);

%%
Y=wm_t;
X = linspace(-1/(2*dkt),1/(2*dkt),length(kt))';

X = 2*pi*X;
dX = X(2)-X(1);


 Yfft = fftshift(ifft(Y));  
        % Remove momentum associated with the sampling frequency
        Yfft = Yfft.*exp(1i*pi/dX*X);    
        % Normalize to one
        Yfft = Yfft/sqrt(sum(Yfft.*conj(Yfft)));   
        
        Yfft = Yfft/sqrt(trapz(X,Yfft.*Yfft));
        
        figure(49);
        clf
        plot(X/pi,Yfft.*X);xlim([-5 5]);
        hold on
        plot(X/pi,sin(X)./X/sqrt(pi));


%%
XL=[.99 1.01];

figure(28);
clf

subplot(121);
plot(k,conj(real(wmk)).*imag(dwnk),'.-')
hold on
plot(kt,conj(real(wm_t)).*imag(dwn_t),'.-')
 xlim(XL)
 
 subplot(122);
plot(k,conj(real(wmk)).*imag(dwnk),'.-')
hold on
plot(kt,conj(real(wm_t)).*imag(dwn_t),'.-')
 xlim(-flip(XL))