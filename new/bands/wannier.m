function [npt] = wannier(npt,opts)
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
% w(x-x0) = integral_k_FBZ exp(i*k*x-x0) * u_{k}(x)*exp(i*theta(k))
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
% Thus finding the wannier function amounts to a inverse Fourier transform 
% of the plan wave cofficients c(l,k).


bands = opts.bands;
% npt.Bands = bands;

%% Extended k-space vector

% Construct extended momentum states
kVec = npt.K; 
kVec = kVec(:);
pvals = diag(npt.PMatrix);
pvals = pvals(:);
 % plane wave index x numK
pall = repmat(pvals,[1 length(kVec)]);
kall = repmat(kVec,[1 length(pvals)])';
keff = pall+kall;
% keff0 = keff;
keff_flat = reshape(keff',[],1);
[keff_flat_ordered,inds_order] = sort(keff_flat);
[keff,inds_unique] = unique(keff_flat_ordered);
% K_eff = keff;

%% Position Vector

% Convert 2*pi*k*x FT convention to k*x
dK = npt.K(2) - npt.K(1);
X = linspace(-1/(2*dK),1/(2*dK),length(keff));

X = 2*pi*X;
dX = X(2)-X(1);
% Assign output X vector
npt.numX = length(X);



%% Calculate Wannier Functions in k-space

for ii = 1:length(npt.depth)    
    depth = npt.depth(ii);
    fprintf(['computing wannier ' ...
        '(U=' num2str(depth) 'Er,' ...
        'Nk = ' num2str(npt.numK) ',' ...
        'Nstates = ' num2str(npt.numStates) ',' ...
        'bands = (' num2str(bands) ...
        ') ...']);    
    tstart = now;
    wannier_k = zeros(length(keff),length(bands));
    wannier_x = zeros(length(keff),length(bands));
    for nn=1:length(bands)
        % Select the band
        n=bands(nn);
        % Get vectors of this band for all quasimomentum
        % l x n x k x ii(plane wave, band, momentum,depth)
        A = npt.bandEigenVectors(:,n,:,ii); 
        % Reorder to a nice 2D matrix
        B = permute(A,[1 3 2]);   % plane wave index x numK
        Bflat = reshape(B',[],1);
        Bordered = Bflat(inds_order);
        Bunique = Bordered(inds_unique);
        Y = Bunique;
        % Normalize to one
        Y = Y/sqrt(sum(Y.*conj(Y)));        
        % Assign k space function
        wannier_k(:,nn)=Y;          
        % Perform FFT to get spatial domain
        Yfft = fftshift(ifft(Y));  
        % Remove momentum associated with the sampling frequency
        Yfft = Yfft.*exp(1i*pi/dX*X');    
        % Normalize to one
        Yfft = Yfft/sqrt(sum(Yfft.*conj(Yfft)));   
        wannier_x(:,nn)=Yfft;
    end
    
    % Assign wannier functions
    npt.K_extended = keff;
    npt.Wannier_K(:,:,ii) = wannier_k;

    npt.X_extended = X/pi;
    npt.Wannier_X(:,:,ii) = wannier_x;     

    %% Harmonic
    % Calculate the harmonic wavefunctions using Hermite Polynomials

    wannier_x_harmonic = zeros(length(X),length(bands));
    for nn=1:length(bands)
        hk = HermitePoly(bands(nn)-1);
        Y2 = X*depth^(1/4);
        psi_ho = polyval(hk,Y2).*exp(-Y2.^2/2);
        psi_ho = psi_ho/sqrt(sum(psi_ho.*conj(psi_ho)));
        wannier_x_harmonic(:,nn) = psi_ho';
    end

    % Assign harmonic wavefunctions
    npt.Wannier_X_Harmonic(:,:,ii) = wannier_x_harmonic;
        
    % Harmonic oscillator length in units of lattice spacing
    npt.Harmonic_Length(ii) = (1/depth)^(1/4)/pi;
    
    tend = toc;
    disp([' done (' num2str(round((tend),3)) ' s)']);    
    

    
    if opts.doPlot
        fprintf(' plotting ... ');
        tstart=now;

        hF_W_probability=figure(10000+ii);
        hF_W_probability.Color='w';
        hF_W_probability.Position=[250 50 800 400];
        hF_W_probability.Name = 'wannier_probability';
        clf

        t=uicontrol('style','text','string',[num2str(depth) 'Er'],...
            'fontname','times','fontsize',12,'backgroundcolor','w');
        t.Position(3:4)=t.Extent(3:4);
        t.Position(1:2)=[5 hF_W_probability.Position(4)-t.Extent(4)];
        for nn=1:length(bands)
            subplot(2,length(bands),nn)
            psi_wannier_k = wannier_k(:,nn);
            rho_wannier_k = psi_wannier_k.*conj(psi_wannier_k);
            rho_wannier_k = rho_wannier_k/max(rho_wannier_k);

            plot(K_eff,rho_wannier_k,'-','linewidth',1);
            xlabel('momentum ($\hbar k_L$)','interpreter','latex');
            str=['$|w_' num2str(bands(nn)) '(k)|^2$'];
            text(.01,.98,str,'units','normalized','interpreter','latex',...
                'verticalalignment','top','fontsize',10);
            set(gca,'xgrid','on','ygrid','on','fontname','times',...     
                'fontsize',8,'yticklabel',{});
            xlim([-20 20]);

            subplot(2,length(bands),length(bands)+nn)
            psi_wannier = wannier_x(:,nn);
            rho_wannier = psi_wannier.*conj(psi_wannier);
            rho_wannier = rho_wannier/max(rho_wannier);

            plot(X/pi,rho_wannier,'-','linewidth',1);
            xlabel('position (site)','interpreter','latex');
            str=['$|w_' num2str(bands(nn)) '(x)|^2$'];

            text(.01,.98,str,'units','normalized','interpreter','latex',...
                'verticalalignment','top','fontsize',10);
            set(gca,'xgrid','on','ygrid','on','fontname','times',...     
                'fontsize',8,'yticklabel',{});
            xlim([-.5 .5]);
            hold on    

            psi_harmonic_x = wannier_x_harmonic(:,nn);
            rho_ho = psi_harmonic_x.*conj(psi_harmonic_x);
            rho_ho = rho_ho/max(rho_ho);
        %     plot(X/pi,rho_ho,'k--','linewidth',1);    
            plot(X/pi,rho_ho,'-','linewidth',2,'color',[.3 .3 .3 .7]);  

        end

        %% Wavefunction Plot
        hF_W_psi=figure(10001);
        hF_W_psi.Color='w';
        hF_W_psi.Position=[1055 50 800 400];
        hF_W_psi.Name = 'wannier_wavefunction';
        clf
        t=uicontrol('style','text','string',[num2str(depth) 'Er'],...
            'fontname','times','fontsize',12,'backgroundcolor','w');
        t.Position(3:4)=t.Extent(3:4);
        t.Position(1:2)=[5 hF_W_psi.Position(4)-t.Extent(4)];
        for nn=1:length(bands)
            subplot(2,length(bands),nn)
            psi_wannier_k = wannier_k(:,nn);

            plot(K_eff,real(psi_wannier_k),'-','linewidth',1);
            hold on
            plot(K_eff,imag(psi_wannier_k),'-','linewidth',1);

            xlabel('momentum ($\hbar k_L$)','interpreter','latex');
            str=['$w_' num2str(bands(nn)) '(k)$'];
            text(.01,.98,str,'units','normalized','interpreter','latex',...
                'verticalalignment','top','fontsize',10);
            set(gca,'xgrid','on','ygrid','on','fontname','times',...     
                'fontsize',8,'yticklabel',{});
            xlim([-20 20]);

            subplot(2,length(bands),length(bands)+nn)
            psi_wannier_x = wannier_x(:,nn);


            plot(X/pi,real(psi_wannier_x),'-','linewidth',1);
            hold on
            plot(X/pi,imag(psi_wannier_x),'-','linewidth',1);

            xlabel('position (site)','interpreter','latex');
            str=['$w_' num2str(bands(nn)) '(x)$'];

            text(.01,.98,str,'units','normalized','interpreter','latex',...
                'verticalalignment','top','fontsize',10);
            set(gca,'xgrid','on','ygrid','on','fontname','times',...     
                'fontsize',8,'yticklabel',{});
            xlim([-.5 .5]);
            hold on    

            psi_harmonic_x = wannier_x_harmonic(:,nn);
            plot(X/pi,psi_harmonic_x,'-','linewidth',2,'color',[.3 .3 .3 .7]);  
        end

        tend = now;
        fprintf([' done (' num2str(round((tend-tstart)*24*60*60,3)) ' s)']);
    else 
        hF_W_probability=[];
        hF_W_psi=[];
    end
    disp(' ');


end

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

function yk = fixGaugeWannier(Cmat,p)
% Cmat is plane wave state (0,1,-1,2,-2,...) x momentum sized (-1,.99,...1)
% yk0 is the input wannier function
% p is the parity
%
% Make it all real for simplicity
%  p = 0,2,4,6 => even parity, real
%  0 ==> int(w)>0
%  2 ==> int(w) <0
%  4 ==> int(w) >0

% p = 1,3,5,7 ==> odd parity, imaginary
%
end




