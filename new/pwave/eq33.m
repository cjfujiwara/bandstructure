% Right hand side of equation 33. Convert energy to d^3/a^3
function d3_a3 = eq33(E,n)

d3_a3 = 0;
for k=0:(n-1)
    d3_a3 = d3_a3 + gamma((k+0.5)/n-0.5*E+0.75)./gamma((k+0.5)/n-E/2-0.75);    
end
d3_a3 = -8./n*d3_a3;

end