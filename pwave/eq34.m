% Right hand side of equation 34. Convert energy to d^3/a^3
function d3_a3 = eq34(E,n)

d3_a3 = 0;
for k=0:(n-1)
    for l=0:(n-1)
        d3_a3 = d3_a3 + ...
            gamma((k+l+1)/n - 0.5*E + 0.25)./gamma((k+l+1)/n - 0.5*E - 1.25);
    end
end
d3_a3 = -8./n^2*d3_a3;

end