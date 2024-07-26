function [npt] = calculateWannierHarmonicCoupling(npt)

dk = npt.K(2)-npt.K(1);

nI = 21;                        % Max site number to couple

iVec = -nI:1:nI;

nB = length(npt.Bands);         % number of bands
nU = size(npt.bandEigenValue,3); % lattice depth vector

curveMat = zeros(nB,nB,2*nI+1,nU);
mat0 = zeros(nB,nB,2*nI+1,nU);

mat1 = zeros(nB,nB,2*nI+1,nU);
mat2 =  zeros(nB,nB,2*nI+1,nU);
mat4 =  zeros(nB,nB,2*nI+1,nU);

K = npt.K_extended;
K=K(:);

% Calculate tunneling elements for each depth
for iU=1:nU
  for iBr = 1:nB
      wm = npt.Wannier_K(:,iBr,iU); % wannier function in momentum

      if mod(iBr-1,2) == 0
        wm = real(wm);
      else
         wm = 1i*imag(wm);
      end
      wm = wm./norm(wm);

      wm = wm(:);
      dwm = gradient(wm,dk);        % gradient of wannier function in momentum    
        
      
      for iBc = 1:nB
            wn = npt.Wannier_K(:,iBc,iU); % wannier function in momentum
             if mod(iBc-1,2) == 0
                wn = real(wn);
              else
                 wn = 1i*imag(wn);
              end
              wn = wn./norm(wn);


            wn = wn(:);
            dwn = gradient(wn,dk);        % gradient of wannier function in momentum
            
            ddwn = gradient(dwn,dk);

            y0 = conj(wm).*wn;
            y1 = conj(wm).*(1i*dwn);

            y2 = wn.*(-1i*conj(dwm));
            y3 = conj(dwm).*dwn;

            y4 = -conj(wm).*ddwn;
            for ii = 1:length(iVec)
                dI = iVec(ii);
                osc = exp(-1i*K*dI*pi);

                v0 = osc.*y0;
                s0 = trapz(v0);
                s0 = real(s0);

                v1 = osc.*y1;
                s1 = trapz(v1);
                s1 = real(s1);

                v2 = osc.*y2;
                s2 = trapz(v2);
                s2 = real(s2);

                v3 = osc.*y3;
                s3 = trapz(v3);
                s3 = real(s3);   

                v4 = osc.*y4;
                s4 = trapz(v4);
                s4 = real(s4);
         
                mat0(iBr,iBc,ii,iU) = s0;
                mat1(iBr,iBc,ii,iU) = s1;
                mat2(iBr,iBc,ii,iU) = s2;

                curveMat(iBr,iBc,ii,iU) = s3;
                mat4(iBr,iBc,ii,iU) = s4;

                % if iBc==2  && dI == 3;
                %     keyboard
                % 
                % end
            end       

      end
  end
end

% npt.WannierHarmonicCoupling = tunneling;
% npt.bandEigenValueAverage = Ebar;
npt.m0 = mat0;

npt.m1 = mat1;
npt.m2 = mat2;
npt.m4 = mat4;

npt.curve = curveMat;


end

