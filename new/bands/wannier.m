function [npt,hF] = wannier(npt)
%% Wannier Function
% The wannier function is a position basis description of states in the
% lattice potential.  They are constructed by an fourier transform of the
% Bloch waves of a particle band over the BZ.
%
% w_n(x)=V/(2*pi)^3 Integral(psi_nk(x),-pi/2,pi/2) dK
%
% We numerically construct this by summing over the plane wave basis
% states.

out=struct;

% Band index
bands=[1 2];

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

for nn=1:length(bands)
    n=bands(nn);

    % Wannier Fucntion is FFT of Bloch waves
    wData=zeros(1,length(X));
    for kk=1:npt.numK
        % Get basis vectors
        v=npt.bandEigenVectors(:,:,kk);

        % Get weights for this band index
        vn=v(:,n);

        % Construct the Bloch wave
        u_nk   = makeU_nk(vn);
        psi_nk = u_nk(X).*exp(1i*K(kk)*X).*exp(-1i * K(kk) * pi*nSite); 
        
        % Normalize to number of sites and for discretization
        psi_nk = psi_nk/(sqrt(Nsites/dL));

        % Normalize the wavefunction
        N0=sum(psi_nk.*conj(psi_nk));
        psi_nk=psi_nk/sqrt(N0);

        % Add this Bloch wave to the Wannier Function
        wData=wData+psi_nk;        
    end

    % Normalize to number of K-states and for discretization
    wData=wData/sqrt(length(K)/dK);

    % Ensure normalization
    N0=sum(wData.*conj(wData));   
    wData=wData/sqrt(N0);    

    switch n
        case 1
            out.WannierS=wData;  
        case 2
            out.WannierP=wData;            
        case 3
            out.WannierD=wData;
        case 4
            out.WannierF=wData;
    end

end

%% Harmonic Oscillator Wavefunctions
Y=X.*depth^(1/4);

HOs=exp(-Y.^2/2);    
Ns=sum(HOs.*conj(HOs));  
HOs=HOs/sqrt(Ns);

HOp=Y.*exp(-Y.^2/2);   
Np=sum(HOp.*conj(HOp));  
HOp=HOp/sqrt(Np);

HOd=(2*Y.^2-1).*exp(-Y.^2/2);  
Nd=sum(HOd.*conj(HOd));  
HOd=HOd/sqrt(Nd);

HOf=(2*Y.^3-3*Y).*exp(-Y.^2/2);    
Nf=sum(HOf.*conj(HOf));  
HOf=HOf/sqrt(Nf);

out.HarmonicS=HOs;
out.HarmonicP=HOp;
out.HarmonicD=HOd;
out.HarmonicF=HOf;
out.X=X;
%% Tunneling
for nn=1:length(bands)
       
    switch bands(nn)
        case 1
            w=out.WannierS;
        case 2
            w=out.WannierP;
        case 3
            w=out.WannierD;
        case 4
            w=out.WannierF;
    end
        
end


npt.Wannier=out;


%% Plot the probability density
hF=figure(223);
set(hF,'color','w','name','wannier');
hF.Position=[320 50 1000 250];
clf


for nn=1:length(bands)
    subplot(1,length(bands),nn);
    
    switch bands(nn)
        case 1
            w=out.WannierS;
            ho=out.HarmonicS;
        case 2
            w=out.WannierP;
            ho=out.HarmonicP;
        case 3
            w=out.WannierD;
            ho=out.HarmonicD;
        case 4
            w=out.WannierF;
            ho=out.HarmonicF;
    end
        
    
    pW=plot(X/pi,real(w.*conj(w)),'LineWidth',2); 
    hold on
    pH=plot(X/pi,real(ho.*conj(ho)),'LineWidth',2); 
    

    xlabel('position (site)');
    ylabel('probability density (arb)');
    legend([pW,pH],{'wannier','harmonic'});

    str=['$' num2str(depth) '~E_R$'];
    text(5,5,str,'interpreter','latex','units','pixels',...
        'verticalalignment','bottom','fontsize',16);
    set(gca,'box','on','linewidth',1,'xgrid','on','ygrid','on');
end


end

