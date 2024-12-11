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
        
for iU=1:nU
    for iBr = 1:nB
        wm = real(npt.Wannier_X(:,iBr,iU));        
        wm = wm/sqrt(trapz(x,wm.*wm));
        
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
        
        
        for iBc = 1:nB
              wn = real(npt.Wannier_X(:,iBc,iU));        
                wn = wn/sqrt(trapz(x,wn.*wn));

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
                
%                 nn=1;
                
                for nn=1:delta    
                   OSC = exp(2*1i*k*nn*pi);
                    D1(iBr,iBc,nn,iU)=trapz(k,C.*OSC);
                end
            % delta=0;
% B=exp(2*1i*k*delta);
%   A=(conj(wkm).*dwkn+conj(dwkm).*wn)*0.5;

            
              if iBc==2 && iBr==1
                   keyboard 
                end
            
            for nn=1:delta
                shift_n = nn*dN;                
                    
                y1 = circshift(wm,shift_n);
                y2 = wn;
                x_shift = x-nn/2;                
                D(iBr,iBc,nn,iU)=trapz(x_shift,y1.*y2.*x_shift);          
                
             
            end
        end
    end
end
      
 

npt.WannierDipoleCoupling = D;

end

