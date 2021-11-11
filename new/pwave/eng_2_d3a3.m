function eng_2_d3a3(n)
%TRAPPED_INTERACTION Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
   n = 5; 
end

%% Equation 33

%% Equation 34


% Ebounds = [-10 -0.225 0.49 ];

Np=1e3;
Ebounds = [-3 1];
a3d3_Mat = zeros(length(Ebounds)-1,Np);
funcs={};
for jj=1:size(Ebounds,1)
    Evec=linspace(Ebounds(jj,1)+1E-5,Ebounds(jj,2)-0.02,Np);

    d3a3 = eq34(Evec,n);
    
   
    
    foo = @(a_in) interp1(d3a3,Evec-1.5,a_in);
    funcs{jj}=foo;    
    a3d3_Mat(jj,:)=d3a3;
end

figure(1);
clf
scatter(1./a3d3_Mat,Evec);
xlim([-100 100]);

Ebounds = [...
    -10 -0.23;
    -0.22 0.574];
a3d3_Mat = zeros(length(Ebounds)-1,Np);
funcs={};
for jj=1:size(Ebounds,1)
    Evec=linspace(Ebounds(jj,1)+1E-5,Ebounds(jj,2),Np);

    d3a3 = eq34(Evec,n);
    
   
    
    foo = @(a_in) interp1(d3a3,Evec-1.5,a_in);
    funcs{jj}=foo;    
    a3d3_Mat(jj,:)=d3a3;
end
% 
% figure(3);
% co=get(gca,'colororder');
% clf
for kk=1:size(a3d3_Mat)
   scatter(1./a3d3_Mat(kk,:), linspace(Ebounds(kk,1)+1E-5,Ebounds(kk,2)-0.02,Np),10,co(kk,:));
   hold on
end
% xlim([-100 100]);

end

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

% Right hand side of equation 33. Convert energy to d^3/a^3
function d3_a3 = eq33(E,n)

d3_a3 = 0;
for k=0:(n-1)
    d3_a3 = d3_a3 + gamma((k+0.5)/n-0.5*E+0.75)./gamma((k+0.5)/n-E/2-0.25);    
end
d3_a3 = -8./n*d3_a3;

end