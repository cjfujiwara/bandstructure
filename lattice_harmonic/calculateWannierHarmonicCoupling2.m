function npt = calculateWannierHarmonicCoupling2(npt)


delta = 51;
nB = length(npt.WannierBands);         % number of bands
nU = size(npt.bandEigenValue,3); % lattice depth vector

x = npt.X_extended';
dx = x(2)-x(1);
dN = round(1/dx);
D = zeros(nB,nB,delta,nU);
% Calculate tunneling elements for each depth

        
k = npt.K_extended;
dk = k(2)-k(1);

fprintf('Calculating dipole coupling ...');
        
% Iterate over all lattice depths
for iU=1:nU
    % Iterate over all rows of bands
    for iBr = 1:nB
        % Momentum Domain
        wmk = npt.Wannier_K(:,iBr,iU);
        if mod(iBr,2)
            wmk = real(wmk);
            wmk = wmk/sqrt(trapz(k,conj(wmk).*wmk));    
            dwmk = gradient(wmk,dk);
        else
            wmk = imag(wmk);
            wmk = wmk/sqrt(trapz(k,conj(wmk).*wmk));    
            dwmk = 1i*gradient(wmk,dk);
            wmk = 1i*wmk;
        end
        % Position Domain
%         wm = real(npt.Wannier_X(:,iBr,iU));        
%         wm = wm/sqrt(trapz(x,wm.*wm));        
        % Iterate over all columns of bands
        for iBc = 1:nB  
            wnk = npt.Wannier_K(:,iBc,iU);
            if mod(iBc,2)
                wnk = real(wnk);
                wnk = wnk/sqrt(trapz(k,conj(wnk).*wnk));    
                dwnk = gradient(wnk,dk);
            else
                wnk = imag(wnk);
                wnk = wnk/sqrt(trapz(k,conj(wnk).*wnk));    
                dwnk = 1i*gradient(wnk,dk);
                wnk = 1i*wnk;
            end
            A = conj(wmk).*dwnk;
            B = conj(dwmk).*wnk;
            C = (A-B)*0.5*1i;
            for nn=1:delta    
                OSC = exp(2*1i*k*(nn*pi/2));
                
          
                D(iBr,iBc,nn,iU)=trapz(k,C.*OSC);
            end
                
            % Position Domain
%              wn = real(npt.Wannier_X(:,iBc,iU));        
%                 wn = wn/sqrt(trapz(x,wn.*wn));            
%             for nn=1:delta
%                 shift_n = nn*dN;         
%                 y1 = circshift(wm,shift_n);
%                 y2 = wn;
%                 x_shift = x-nn/2;                
%                 D(iBr,iBc,nn,iU)=trapz(x_shift,y1.*y2.*x_shift);   
%             end
        end
    end
end   

disp('done ...');

 npt.WannierDipoleCoupling = D;
end

