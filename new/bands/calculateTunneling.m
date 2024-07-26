function [npt,hF] = calculateTunneling(npt)
    
dk = npt.K(2)-npt.K(1);
deltaImax = 101;
tunneling = zeros(size(npt.bandEigenValue,1),deltaImax,size(npt.bandEigenValue,3));
Ebar = zeros(size(npt.bandEigenValue,1),size(npt.bandEigenValue,3));
% Calculate tunneling elements for each depth
for ind_depth=1:size(npt.bandEigenValue,3)
  for ind_band = 1:size(npt.bandEigenValue,1)
        Eq = npt.bandEigenValue(ind_band,:,ind_depth);
        Eq = Eq(:);

        Ebar(ind_band,ind_depth) = 0.5*trapz(Eq)*dk; % average band energy

        for ii = 1:deltaImax
            yy = cos(pi*npt.K(:)*ii);
            t = -0.5*trapz(yy.*Eq)*dk;
            tunneling(ind_band,ii,ind_depth) = t;     
        end         
  end
end

npt.Tunneling = tunneling;
npt.bandEigenValueAverage = Ebar;



%% Band Curvature

% Band curavture at Gamma point (center of FBZ)
npt.BandCurvatureG = zeros(size(npt.bandEigenValue,1),size(npt.bandEigenValue,3));
% Band curvature at X point (edge of FBZ)
npt.BandCurvatureX = zeros(size(npt.bandEigenValue,1),size(npt.bandEigenValue,3));
for ind_depth=1:size(npt.bandEigenValue,3)
  for ind_band = 1:size(npt.bandEigenValue,1)
      t_band = tunneling(ind_band,:,ind_depth);
      t_band = t_band(:);
      n_vec = 1:deltaImax;
      n_vec =n_vec(:);
      kappa = pi^2*sum(n_vec.^2.*t_band);
      kappa = min([1 kappa]); % tunneling expression fails at very low lattice depths, cap it
      npt.BandCurvatureFromTunneling(ind_band,ind_depth)=kappa;

      uu = npt.bandEigenValue(ind_band,:);uu=uu(:);
      kk = npt.K;kk=kk(:);

      [~,ind]=min(abs(kk));
      kL = kk(ind-1);uL =uu(ind-1);
      kC = kk(ind);uC =uu(ind);
      kR = kk(ind+1);uR =uu(ind+1);
      pp=polyfit([kL kC kR],[uL uC uR],2);
      npt.BandCurvatureG(ind_band,ind_depth)=pp(1);

      % [~,ind]=max(kk);
      
        k1 = kk(1);
        k2 = kk(5);
        u1 = uu(1);
        u2 = uu(5);
        kappa = (u2-u1)/(k2-k1)^2;
      % kL = kk(ind-2);uL =uu(ind-2);
      % kC = kk(ind-1);uC =uu(ind-1);
      % kR = kk(ind);uR =uu(ind);
      % myfit = 
      pp=polyfit(kk(1:10),uu(1:10),2);

      npt.BandCurvatureX(ind_band,ind_depth)=kappa;
  end
end

end

