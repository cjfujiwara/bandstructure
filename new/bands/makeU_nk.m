function [u_nk] = makeU_nk(v)
% Given an input vector v of plane wave coefficients, return the function 
% u_{n,k}(x).

    % 1/sqrt(pi), exp(2*x*1i)/sqrt(pi),exp(-2*x*1i)/sqrt(pi),
    % exp(4*x*1i)/sqrt(pi),exp(-4*x*1i)/sqrt(pi) ...

    % |0>
    u_nk=@(phi) v(1)/sqrt(pi);    
    nMax=(length(v)-1)/2;
    
    % |n>
    for n=1:nMax
       u_nk=@(phi) u_nk(phi) + ...
           v(2*n+0)*exp(+2*1i*n*phi)/sqrt(pi)+...
           v(2*n+1)*exp(-2*1i*n*phi)/sqrt(pi);    
    end
end
