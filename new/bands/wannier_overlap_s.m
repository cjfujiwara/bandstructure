function npt = wannier_overlap_s(npt)
%WANNIER_OVERLAP Summary of this function goes here
%   Detailed explanation goes here
%
%
% Some useful transformations
%
% a_L      = pi/k_L = lambda/2
% l_ho/a_L = 1/(k_L*depth^(1/4))
%          = (1/depth)^(1/4)*1/pi
%
% To first order the interaction energy for two fermions in the lowest band
% interacting with an s-wave scattering length a is given by:
%
% U := g Integral |w(x)|^4 d^3x = g \prod_i |w(x_i)|^4 dx_i
%
% where g = 4*pi*a*hbar^2/m
%
% It can be shown that
%
% g/Er = (a/l_ho) * depth^(-1/4) & (8*pi/k_L^3)
%
% This can be expressed in terms of recoil energy and the dimensionless
% position X:=k_L x
%
% Integral |w(x)|^4 dx = k_L Integral |w(X)|^4 dX
%     
% Subject to unitarity |<w(x)|w(x)>|^2=|<w(X)|w(X)>|^2 = 1
%
% So that the interaction energy may be written as:
%
% U/Er : = g/Er k_L^3 (Integral |w(X)|^4 dX)^3
%      : = (a/l_ho) * depth^(-1/4) * 8*pi * (Integral |w(X)|^4 dX)^3
%        = (a/l_ho) * depth^(-1/4) * 8*pi * Is^3
%
% Assuming a harmonic wavefunction
%
% Is = depth^(1/4)/sqrt(2*pi)
%
% So that the interaction energy simplifies to
%
% U/Er : = (a/l_ho) sqrt(2/pi) * sqrt(4*u)
% U/Er : = (a/l_ho) sqrt(2/pi) * (hbar*omega/Er)
% 
% And we recover the linear expectation
%
% U = (a/l_ho) sqrt(2/pi) * (hbar*omega)

disp('Calculating wannier function overlap integrals');

%% Check the Harmonic Integrals

OverlapS = struct;
OverlapS.depth = npt.depth;

Is_ho_mat = zeros(length(npt.Bands),length(npt.depth));
Is_ho_th_mat = zeros(length(npt.Bands),length(npt.depth));
Is_wn_mat = zeros(length(npt.Bands),length(npt.depth));

U_ho_mat = zeros(length(npt.Bands),length(npt.depth));
U_ho_th_mat = zeros(length(npt.Bands),length(npt.depth));
U_wn_mat = zeros(length(npt.Bands),length(npt.depth));

disp('Checking s-band overlap integrals');
for ii=1:length(npt.depth)
    depth = npt.depth(ii);
    for nn=1:length(npt.Bands)
        % Grab the position vector
        x = npt.X;dx = x(2) - x(1);

        % Harmonic oscillattor length scale (to truncate integration)
        l_ho = pi*npt.Harmonic_Length(ii);   

        % Grab the wavefunctions and scale
        y_ho = npt.Wannier_X_Harmonic(:,nn,ii)/sqrt(dx);
        y_wn = npt.Wannier_X(:,nn,ii)/sqrt(dx);

        % Sanity check normalization
        %        fun = @(x_in) interp1(x,y.^2,x_in);

        % Interaction
        fun_ho = @(x_in) interp1(x,(y_ho.*conj(y_ho)).^2,x_in);
        fun_wn = @(x_in) interp1(x,(y_wn.*conj(y_wn)).^2,x_in);

        % Calculate integral
        Is_ho = integral(fun_ho,-10*l_ho,10*l_ho);
        Is_wn = integral(fun_wn,-10*l_ho,10*l_ho);

        Is_ho_th = depth^(1/4)/sqrt(2*pi);

        Is_ho_mat(nn,ii) = Is_ho;
        Is_ho_th_mat(nn,ii) = Is_ho_th;
        Is_wn_mat(nn,ii) = Is_wn;

        U_ho_mat(nn,ii) = depth^(-1/4) * 8*pi * Is_ho^3;
        U_ho_th_mat(nn,ii) = depth^(-1/4) * 8*pi * Is_ho_th^3;
        U_wn_mat(nn,ii) = depth^(-1/4) * 8*pi * Is_wn^3;

        str = [num2str(depth) ' Er : ' ...
            'HO Theory = ' num2str(round(Is_ho_th,3)) ',' ...
            'HO Numeric = ' num2str(round(Is_ho,3)) ',' ...
            'Wannier = ' num2str(round(Is_wn,3))];
        disp(str);
    end
end

OverlapS.Is_Harmonic = Is_ho_mat;
OverlapS.Is_Wannier = Is_wn_mat;
OverlapS.Is_Harmonic_Theory = Is_ho_th_mat;

OverlapS.U_Harmonic_Er = U_ho_mat;
OverlapS.U_Wannier_Er = U_wn_mat;
OverlapS.U_Harmonic_Theory_Er = U_ho_th_mat;


npt.OverlapS = OverlapS;
end

