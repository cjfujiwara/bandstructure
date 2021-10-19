function [npt,hF_W_k] = wannier(npt)
%% Wannier Function
% The wannier function is a position basis description of states in the
% lattice potential.  They are constructed by an fourier transform of the
% Bloch waves of a particle band over the BZ.
%
% https://en.wikipedia.org/wiki/Wannier_function
%
% w_n(x)=V/(2*pi)^3 Integral(psi_nk(x),-pi/2,pi/2) dK
%
% We numerically construct this by summing over the plane wave basis
% states.
%
% This can be simplfied because our basis is already the fourier basis.
% Drop the band label for notational simplicity
%
% w(x-x0) = integral_k_FBZ exp(i*k*x-x0) * u_{k}(x)
%           = integral_k_FBZ exp(i*k*x-x0) * \sum_l c(l,k) exp(2*i*l*x)
%           = integral_k_FBZ \sum_l c(l,k) exp(i*(k+2*l)*x)*exp(-i*k*x0)
%
% Here note that the coefficients for the nth band eigenstates c(l,k) are
% only defined within the FBZ and exp(-i*k*x0) is periodic with respect to 
% a FBZ.  We can simplify by noting that 2*l+k spans the entire k space 
% uniquely so that the integral simplifies to
%
% w(x-x0) = integral_K exp(i*K*x-x0) c_eff(K)
%
% Where c_eff(K) is the unfolded c(l,k) where K spans the entirety of
% momentum space so that k = mod(K,FBZ), K = k+2*l for a given BZ l.
%
% Thus finding the wannier function amounts to a Fourier transform of the
% cofficients c(l,k) in the plane wave basis.

out=struct;

% Band index
bands=[1 2 3 4];

% Central site
nSite=0;

% Quasimomentum vector
K=npt.K;
dK=K(2)-K(1);

% Position vector
X=npt.X;
dX=X(2)-X(1);

% Lattice depth
depth=npt.depth;

% Number of sites
Nsites=(X(end)-X(1))/pi;

% Point spacing
dL=X(2)-X(1);

%% Calculate wannier functions
% 
% tic
% for nn=1:length(bands)
%     n=bands(nn);
% 
%     % Wannier Fucntion is FFT of Bloch waves
%     wData=zeros(1,length(X));
%     for kk=1:npt.numK
%         % Get basis vectors
%         v=npt.bandEigenVectors(:,:,kk);
% 
%         % Get weights for this band index
%         vn=v(:,n);      
%         
%         % Construct the Bloch wave
%         u_nk   = makeU_nk(vn);
%         psi_nk = u_nk(X).*exp(1i*K(kk)*X).*exp(-1i * K(kk) * pi*nSite); 
%         
%         % Normalize to number of sites and for discretization
%         psi_nk = psi_nk/(sqrt(Nsites/dL));
% 
%         % Normalize the wavefunction
%         N0=sum(psi_nk.*conj(psi_nk));
%         psi_nk=psi_nk/sqrt(N0);
% 
%         % Add this Bloch wave to the Wannier Function
%         wData=wData+psi_nk;        
%     end
% 
%     % Normalize to number of K-states and for discretization
%     wData=wData/sqrt(length(K)/dK);
% 
%     % Ensure normalization
%     N0=sum(wData.*conj(wData));   
%     wData=wData/sqrt(N0);    
% 
%     switch n
%         case 1
%             out.WannierS=wData;  
%         case 2
%             out.WannierP=wData;            
%         case 3
%             out.WannierD=wData;
%         case 4
%             out.WannierF=wData;
%     end
% 
% end
% 
% toc
%% Calculate Wannier Functions in k-space
fprintf('Calculating wannier functions via fft ...');
tstart = now;

% Find ordering of plane wave basis
v = -(npt.numStates-1)/2:1:(npt.numStates-1)/2;
ind = 2*abs(v)+(v<=0);

% Calculate unfolded brillouin zone
K_eff = repmat(npt.K,[npt.numStates 1]);
K_eff = K_eff + 2*v';K_eff = K_eff';K_eff = K_eff(:);


% Remove duplicates at corners of zones
for jj=1:(npt.numStates-1)
    i1 = (npt.numK-1)*jj+1;
    i2 = i1 + 1;
    K_eff(i2)=[];
end

X_eff = linspace(-1/dK,1/dK,length(K_eff));


wannier_k ={};
for nn=1:length(bands)
    % Select the band
    n=bands(nn);
    
    % Get vectors of this band for all quasimomentum
    % l x n x k (plane wave, band, momentum)
    A = npt.bandEigenVectors(:,n,:); 
    
    % Reorder to a nice 2D matrix
    B = permute(A,[1 3 2]);   
    
    % Rearrange cofficients in lowest to high K_eff
    C = B(ind,:);
    
    % Transform and get values
    D = C';
    Y = D(:);
    % Remove duplicates at corners of zones
    for jj=1:(npt.numStates-1)
        i1 = (npt.numK-1)*jj+1;
        i2 = i1 + 1;        
        
        Y(i1)=(Y(i1)+Y(i2))/2;
        Y(i2)=[];
    end

    % Assign k space function
    wannier_k{nn}=Y;
    
    % Normalize to one
    wannier_k{nn}=wannier_k{nn}/sqrt(sum(wannier_k{nn}.*conj(wannier_k{nn})));
    
    % Perform FFT to get spatial domain
    wannier_x{nn}=fftshift(fft(Y));
    
    % Normalize to one
    wannier_x{nn}=wannier_x{nn}/sqrt(sum(wannier_x{nn}.*conj(wannier_x{nn})));
end

tend = now;
disp(['done (' num2str(round((tend-tstart)*24*60*60,3)) ' s)']);

%% Coefficient plot

hF_W_k=figure(10000);
hF_W_k.Color='w';
hF_W_k.Position=[250 50 900 400];

t=uicontrol('style','text','string',[num2str(depth) 'Er'],...
    'fontname','times','fontsize',12,'backgroundcolor','w');
t.Position(3:4)=t.Extent(3:4);
t.Position(1:2)=[5 hF_W_k.Position(4)-t.Extent(4)];
for nn=1:length(bands)
    subplot(2,length(bands),nn)
    psi_wannier_k = wannier_k{nn};
    rho_wannier_k = psi_wannier_k.*conj(psi_wannier_k);
    rho_wannier_k = rho_wannier_k/max(rho_wannier_k);
    
    plot(K_eff,rho_wannier_k,'-','linewidth',2);
    xlabel('momentum ($\hbar k_L$)','interpreter','latex');
%     ylabel('$|w_n(k)|^2$','interpreter','latex');
    str=['$|w_' num2str(bands(nn)) '(k)|^2$'];
    text(.01,.98,str,'units','normalized','interpreter','latex',...
        'verticalalignment','top','fontsize',10);
    set(gca,'xgrid','on','ygrid','on','fontname','times',...     
        'fontsize',8,'yticklabel',{});
    xlim([-50 50]);
    
    subplot(2,length(bands),length(bands)+nn)
    psi_wannier = wannier_x{nn};
    rho_wannier = psi_wannier.*conj(psi_wannier);
    rho_wannier = rho_wannier/max(rho_wannier);
    
    plot(X_eff,rho_wannier,'-','linewidth',2);
    xlabel('position (site)','interpreter','latex');
%     ylabel('$|w_n(x)|^2$','interpreter','latex');
    str=['$|w_' num2str(bands(nn)) '(x)|^2$'];

    text(.01,.98,str,'units','normalized','interpreter','latex',...
        'verticalalignment','top','fontsize',10);
    set(gca,'xgrid','on','ygrid','on','fontname','times',...     
        'fontsize',8,'yticklabel',{});
    xlim([-.5 .5]);
    hold on
    
    
    hk = HermitePoly(bands(nn)-1);
    X2 = linspace(-5,5,1000);
    Y2 = X2*depth^(1/4);
    psi_ho = polyval(hk,Y2).*exp(-Y2.^2/2);
    rho_ho = psi_ho.*conj(psi_ho);
    rho_ho = rho_ho/max(rho_ho);
    plot(X2/pi,rho_ho,':','linewidth',2);    
%     
%         legend({'w','ho'},'location','southeast');

end



%% Harmonic Oscillator Wavefunctions
% 
% Y=X.*depth^(1/4);
% 
% HOs=exp(-Y.^2/2);    
% Ns=sum(HOs.*conj(HOs));  
% HOs=HOs/sqrt(Ns);
% 
% HOp=Y.*exp(-Y.^2/2);   
% Np=sum(HOp.*conj(HOp));  
% HOp=HOp/sqrt(Np);
% 
% HOd=(2*Y.^2-1).*exp(-Y.^2/2);  
% Nd=sum(HOd.*conj(HOd));  
% HOd=HOd/sqrt(Nd);
% 
% HOf=(2*Y.^3-3*Y).*exp(-Y.^2/2);    
% Nf=sum(HOf.*conj(HOf));  
% HOf=HOf/sqrt(Nf);
% 
% out.HarmonicS=HOs;
% out.HarmonicP=HOp;
% out.HarmonicD=HOd;
% out.HarmonicF=HOf;
% out.X=X;
%% Tunneling
% for nn=1:length(bands)
%        
%     switch bands(nn)
%         case 1
%             w=out.WannierS;
%         case 2
%             w=out.WannierP;
%         case 3
%             w=out.WannierD;
%         case 4
%             w=out.WannierF;
%     end
%         
% end
% 
% 
% npt.Wannier=out;


%% Plot the probability density
% hF=figure(223);
% set(hF,'color','w','name','wannier');
% hF.Position=[250 250 900 200];
% clf
% 
% 
% for nn=1:length(bands)
%     subplot(1,length(bands),nn);
%     
%     switch bands(nn)
%         case 1
%             w=out.WannierS;
%             ho=out.HarmonicS;
%         case 2
%             w=out.WannierP;
%             ho=out.HarmonicP;
%         case 3
%             w=out.WannierD;
%             ho=out.HarmonicD;
%         case 4
%             w=out.WannierF;
%             ho=out.HarmonicF;
%     end
%         
%     
%     pW=plot(X/pi,real(w.*conj(w)),'LineWidth',2); 
%     hold on
%     pH=plot(X/pi,real(ho.*conj(ho)),'LineWidth',2); 
%     
% 
%     xlabel('position (site)');
%     ylabel('probability density (arb)');
%     legend([pW,pH],{'w','HO'});
% 
%     str=['$' num2str(depth) '~E_R$'];
%     text(5,5,str,'interpreter','latex','units','pixels',...
%         'verticalalignment','bottom','fontsize',16);
%     set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
% end


end

function hk = HermitePoly(n)
if n==0 
    hk = 1;
elseif n==1
    hk = [2 0];
else
    
    hkm2 = zeros(1,n+1);
    hkm2(n+1) = 1;
    hkm1 = zeros(1,n+1);
    hkm1(n) = 2;
    for k=2:n
        
        hk = zeros(1,n+1);
        for e=n-k+1:2:n
            hk(e) = 2*(hkm1(e+1) - (k-1)*hkm2(e));
        end
        
        hk(n+1) = -2*(k-1)*hkm2(n+1);
        
        if k<n
            hkm2 = hkm1;
            hkm1 = hk;
        end
        
    end
    
end
end

