function [h_mat,p_mat]=makeHmatrix(nfo)

% Parameters of hamiltonian
if nargin==0
    depth=10;
    k=0;
    numStates=21;    
else
    depth=      nfo.depth;                  % static lattice depth
    k=          nfo.k;                      % quasimomentum
    numStates=  nfo.numStates;              % number of basis states
end


%% Theoretical Description
% We desire to solve the eigenvalue and eigenvector problem of a cosine
% periodic lattice.
%
% This is given by
%
% H = p^2/2m + U(x)
% H = p^2/2m - v0 cos(kL x)^2
%
% Using the unitless operators phi:=kx, p_phi:=p/(hhbar*kL) his reduces to
%
% h = p_phi^2 - v0/2(cos(2 phi)+1)
%
% After applying Bloch's theorem, the eigenstates are deteremined by the
% following form on the periodic domain  phi = [-pi/2, pi/2]
%
% psi = u_{n,k}exp(i k phi)
%
% Which results in solving for the u_{n,k}
%
% e_{n,k} u_{n,k} = (p+k)^2+U(x) u_{n,k}
%
% where k is the quasimomentum [-1,1] and p is dimesional momentum operator
% which in position basis is
% 
% p_phi = -i d/dphi

% The basis are pi periodic plane wave state 
%
%   |0> := 1/sqrt(pi)
%   |n> := exp(2n i phi)/sqrt(pi)
%
% where n = +- {1,2,3,4,...}
%
% Momentum operator (diagonal)
% p|n> = -i d/dphi exp(2n i phi)/sqrt(pi)
% p|n> = 2n exp(2n i phi)/sqrt(pi)
% p|n> = 2n |n>
%
% Cosine operator (off diagonal)
% cos(2 phi)    = 0.5*(exp(+2 i phi) + exp(-2 i phi))
% cos(2 phi)|n> = 0.5*(exp(+2 i phi) + exp(-2 i phi)) exp(2 n i phi)/sqrt(pi)
% cos(2 phi)|n> = 0.5 (exp(2 (n+1) i phi) + exp(2 (n-1) i phi))/sqrt(pi)
% cos(2 phi)|n> = 0.5 (|n+1> + |n-1>)
%
% Potential operator (offdiagonal)
% - v0/2(cos(2 phi)+1) |n> = -v0/4(|n+1>+|n-1>) - v0/2 |n>
%                      
% Basis ordering 
% To maintain simple matrices, we shall order the states as
% {|0>,|1>,|-1>,|2>,|-2>,...,|n>,|-n>}


%% Matrix Construction

nMax=(numStates-1)/2;

% identity matrix
I_mat = eye(numStates);

% momentum matrix/operator
nV = zeros(numStates,1);
for ii=2:2:(numStates-1)
    nV(ii)   = +ii;
    nV(ii+1) = -ii;
end
p_mat = diag(nV);

U_mat = zeros(numStates,numStates);
% potential energy matrix
for n=-nMax:1:nMax
   
    % |0>
    if n==0
        ind=1;
        U_mat(ind,ind)=1/2;
        U_mat(ind+1,ind)=1/4;
        U_mat(ind+2,ind)=1/4;        
    end    
    
    % |1>
    if n==1
        ind=2;
        U_mat(ind,ind)=1/2;     % 1
        U_mat(1,ind)=1/4;       % 0
        U_mat(ind+2,ind)=1/4;   % 2
    end
    
    % |-1>
    if n==-1
        ind=3;
        U_mat(ind,ind)=1/2;   % -1
        U_mat(1,ind)=1/4;     % 0
        U_mat(ind+2,ind)=1/4; % -2
    end
    
    % |n>
    if abs(n)>1 && n>0
        ind=2*n;
        U_mat(ind,ind)=1/2;     % 1
        U_mat(ind-2,ind)=1/4;       % 0
        U_mat(ind+2,ind)=1/4;   % 2
    end
    
    % |-n>
    if abs(n)>1  && n<0       
        ind=2*abs(n)+1;
        U_mat(ind,ind)=1/2;     % 1
        U_mat(ind-2,ind)=1/4;       % 0
        U_mat(ind+2,ind)=1/4;   % 2
    end             
end


% Truncate
U_mat = U_mat(1:numStates,1:numStates);

U_mat(numStates-1,numStates-1) = 3/4;
U_mat(numStates,numStates) = 3/4;

% Form the hamiltonian
h_mat = (p_mat + k*I_mat)^2 - U_mat*depth;

% keyboard

% h_mat = (p_mat + abs(k)*I_mat)^2 - U_mat*depth;

end

