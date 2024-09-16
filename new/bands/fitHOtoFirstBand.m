function [ho_input] = fitHOtoFirstBand(ho_input)

e1 = ho_input.BandRanges(1,:);

i1 = find(ho_input.EigenValues>e1(2),1);
i1 = i1-1;

eng = ho_input.EigenValues(1:i1);
n = 0:(i1-1);

p1=polyfit(n,eng,1);

ho_input.FirstBandLinearSlope = p1(1);
end

