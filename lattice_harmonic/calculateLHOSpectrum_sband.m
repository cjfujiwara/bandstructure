function [npt,output] = calculateLHOSpectrum_sband(npt,opts)
% Calculate the 1D eigenspectrum.  Include tunneling of arbitrary order,
% but only include single band physics

nBands                              = 1;
uu                                  = 1;
n                                   = opts.NumSites;
jjMax                               = opts.MaxTunnelingOrder;
Tmat                                = makeTmatrix(n,jjMax); % Could be made faster
Omega                               = opts.Omega;
fr                                  = npt.fr;

%% Construct change of basis matrix
% Unitary operator to convert form single site states to odd/even pairs

nc                  = (n+1)/2; % Index which is the center
u1                  = eye(n);
u2                  = flip(u1,1);
u1(nc:end,nc:end)   = -u1(nc:end,nc:end);

% Change of basis matrix
U                   = u1+u2;
U                   = U/sqrt(2);  % To preserve the norm
U(nc,nc)            = 1; % Center site is unchanged

% 1:nc are even states
% nc+1:end are odd states

%% Construct Harmonic Position Operator
r = (n-1)/2;        % [-r,r] is the position in size
x = -r:1:r;         % position operator is symmetric
x2 = x.^2;          % position squared

%% Construct Potential Energy Operator
V = Omega*diag(x2);

%% Construct Kinetic Energy Operator

t = npt.Tunneling(:,:,uu)*npt.fr; % Matrix of tunneling elements (band,site)
nn=1; % Band index
% Could be made faster
% Kinetic Energy operator for each band
T = zeros(n,n,nBands);  
for jj = 1:jjMax % Iterate over all tunneling order
    T(:,:,nn) = T(:,:,nn) + t(nn,jj)*Tmat(:,:,jj);           
end

%% Construct Hamiltonian and Separate into Odd/Even

Ebar        = npt.bandEigenValueAverage(nn,uu)*fr;      
E0          =  eye(n)*Ebar;

H           = T + V;                        % s-band Hamiltonian
Hoddeven    = U*H*ctranspose(U);            % odd/even Hamiltonian
Hoddeven(abs(Hoddeven)<1e-3)=0;             % Get rid of numerical noise

% Add band offset
Hoddeven = Hoddeven+E0;
H = H + E0;

% Divide into Odd and Even Hamiltonians
Heven = Hoddeven(1:nc,1:nc);             % even sector
Hodd  = Hoddeven((nc+1):end,(nc+1):end); % odd sector

%% Solve Modified Hamiltonians
% Use these to compare and make sure the odd/even thing works if desired

% Solve Total Hamiltonian (comment for speed)
% [c,eng] = eig(H);eng=diag(eng);

% Solve Odd/Even Hamiltonian (comment for speed)
% [c_oddeven,eng_oddeven]=eig(Hoddeven);eng_oddeven=diag(eng_oddeven);

%% Solve

% Solve Even Hamiltonian
[c_even,eng_even]=eig(Heven);eng_even=diag(eng_even);
s=sign(c_even(end,:));
S = repmat(s,[size(c_even,1) ,1]);
c_even = S.*c_even;

% Solve Odd Hamiltonian
[c_odd,eng_odd]=eig(Hodd);eng_odd=diag(eng_odd);
s=sign(c_odd(1,:));
S = repmat(s,[size(c_odd,1) ,1]);
c_odd = S.*c_odd;


% Interleave even and odd eigenvalues
eng_odd(end+1)=NaN; % Add extra value so same length
eng_oddeven2 = reshape([eng_even';eng_odd'],1,[]);
eng_oddeven2(end)=[]; % remove nan guy
eng_oddeven2=eng_oddeven2';
eng_odd(end)=[];

% Interleave even and odd eigenvectors
c_even_full = [c_even; zeros(size(c_odd,1),size(c_even,2))];
c_odd_full = [zeros(size(c_even,1),size(c_odd,2)); c_odd];   % blk diagonal
c_odd_full(:,end+1)=NaN(n,1);
c_oddeven2 = reshape([c_even_full;c_odd_full],n,[]);
c_oddeven2(:,end)=[];

% Convert into original position basis
c_oddeven2 = U*c_oddeven2;

%% Dipole Moment Operator

xT  = x';
M1  = c_oddeven2;
M2  = M1;
M3  = repmat(xT,[1 size(M2,2)]);
M4  = M2.*M3;
D   = ctranspose(M1)*M4;


%% Initialize ouput

output = struct;
output.Depth = npt.depth;
output.Tunneling = npt.Tunneling;
output.NumSites = n;
output.NumBands = nBands;
output.MaxTunnelingOrder = jjMax;
output.Omega = Omega;   
output.omega = opts.omega;
output.fr = fr;
output.PositionVector = repmat(x,[1 nBands]);
output.EigenValues(:,uu) = eng_oddeven2;
output.EigenVectors(:,:,uu) = c_oddeven2;    
output.DipoleOperator = D;

% output.EigenValues = zeros(output.NumSites*output.NumBands,length(npt.depth));
% output.EigenVectors = zeros(output.NumSites*output.NumBands,output.NumSites*output.NumBands,length(npt.depth));

output.BandRanges = zeros(output.NumBands,2,length(npt.depth));
output.BandProjection = zeros(output.NumSites*output.NumBands,output.NumBands,length(npt.depth));



output.BandProjection(:,nn,uu) = ones(size(output.BandProjection,1),1);
output.BandRanges(nn,1,uu) = min(npt.bandEigenValue(nn,:,uu))*fr; 
output.BandRanges(nn,2,uu) = max(npt.bandEigenValue(nn,:,uu))*fr; 



end
