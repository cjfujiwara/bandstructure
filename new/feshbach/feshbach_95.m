function a= feshbach_95(B)
%FESHBACH_79 Summary of this function goes here
%   Detailed explanation goes here

% PhysRevLett.90.230404
% a_bg = 174;a_bg=167;
% Delta = 9.7;
% B0 = 224.21;

% PhysRevLett.90.230404
% a_bg = 174;a_bg=167;
% Delta = 9.7;
% B0 = 224.21;


% Ludewig CC
Delta = 7.2;%Delta=7.6;
% Delta = 9.7;
B0 = 224.2;
a_bg = 167.3;

a = a_bg * (1 - Delta./(B-B0));
end

