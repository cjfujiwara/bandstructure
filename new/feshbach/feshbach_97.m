function a= feshbach_97(B)
%FESHBACH_79 Summary of this function goes here
%   Detailed explanation goes here
% our paper arxiv 2101.02700
a_bg = 166.978;
Delta = 6.910;
B0 = 202.15;

a = a_bg * (1 - Delta./(B-B0));
end

