
U=linspace(0,200,1e3);

omega = 2*pi*npt.fr*sqrt(4*U);
a_ho = sqrt(npt.hbar./(npt.m.*omega))./npt.a0;
figure(8)
pLinearHO=plot(U,omega/(2*pi) *sqrt(2/pi).*(166.978./a_ho),'k-','linewidth',2);
