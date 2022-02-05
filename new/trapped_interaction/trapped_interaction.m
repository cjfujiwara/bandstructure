function [funcs] = trapped_interaction
%TRAPPED_INTERACTION Summary of this function goes here
%   Detailed explanation goes here

a = @(E) (sqrt(2)*gamma(-E/2+3/4)./gamma(-E/2+1/4)).^(-1);

Ebounds = [-20 0.5:2:10.5];
Np=5e5;

aMat = zeros(length(Ebounds)-1,Np);
funcs={};
for jj=1:length(Ebounds)-1
    Evec=linspace(Ebounds(jj)+1E-5,Ebounds(jj+1)-1E-5,Np);
    a_out = a(Evec);
    foo = @(a_in) interp1(a_out,Evec-1.5,a_in);
    funcs{jj}=foo;    
    aMat(jj,:)=a_out;
end

end

