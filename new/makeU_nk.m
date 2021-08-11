function [u_nk] = makeU_nk(v)
%MAKEWAVEFUNCTION Summary of this function goes here
%   Detailed explanation goes here

    % |0>
    u_nk=@(phi) v(1)*sqrt(1/pi);
    
    nMax=(length(v)-1)/2;
    % |n>
    for n=1:nMax
       u_nk=@(phi) u_nk(phi) + ...
           v(2*n+0)*exp(+2*1i*n*phi)+...
           v(2*n+1)*exp(-2*1i*n*phi);    
    end
end
