function [V_0,V_1] = feshbach_77(B)
%FESHBACH_79 Summary of this function goes here
%   Detailed explanation goes here
% our paper arxiv 2101.02700

% mL = 0
B0_0 = 198.803;
Vbg_0 = -(108.0)^3;
Delta_0 = -19.89;
R0_0 = 49.4;

% mL = 1
B0_1 = 198.300;
Vbg_1 = -(107.35)^3;
Delta_1 = -19.54;
R0_1 = 48.9;

V_0 = Vbg_0 * (1 - Delta_0./(B-B0_0));
V_1 = Vbg_1 * (1 - Delta_1./(B-B0_1));

end

