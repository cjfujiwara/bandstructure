function funcs = eng_2_d3a3
%TRAPPED_INTERACTION Summary of this function goes here
%   Detailed explanation goes here
n=1;
%% Equation 33

% n = 1 spacing is 2 always
% n = 2 spacing is 1.22 first, then to 1
% n = 3 spacing is 1.150
% n = 4
%% Equation 34

Np=1e5;
Ebounds = [-1 2];
a3d3_Mat = zeros(length(Ebounds)-1,Np);
funcs={};
for jj=1:size(Ebounds,1)
    Evec=linspace(Ebounds(jj,1),Ebounds(jj,2),Np);
    d3a3 = eq34(Evec,n);      
    foo = @(a_in) interp1(d3a3,Evec-1.5,a_in);
    funcs{jj}=foo;    
    a3d3_Mat(jj,:)=d3a3;
end

% figure(1);
% clf
% scatter(1./a3d3_Mat,Evec);
% xlim([-100 100]);

Ebounds = [...
    -5 -0.5;
    -0.5 1.5;
    1.5 3.5;
    3.5 5.5;
    5.5 7.5;
    7.5 9.5];

a3d3_Mat = zeros(length(Ebounds)-1,Np);
funcs={};
for jj=1:size(Ebounds,1)
    Evec=linspace(Ebounds(jj,1)+1E-5,Ebounds(jj,2)-1e-5,Np);
    d3a3 = eq34(Evec,n);      
    foo = @(a_in) interp1(1./d3a3,Evec,a_in);
    funcs{jj}=foo;    
    a3d3_Mat(jj,:)=d3a3;
end

% figure(3);
% co=get(gca,'colororder');
% clf
% for kk=1:size(a3d3_Mat)
%    scatter(1./a3d3_Mat(kk,:), linspace(Ebounds(kk,1)+1E-5,Ebounds(kk,2)-0.02,Np),10,co(kk,:));
%    hold on
% end
% xlim([-100 100]);

end

