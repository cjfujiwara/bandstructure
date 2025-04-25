function z = conductivity2(f,T,G,lho)
% Author : CJ Fujiwara
%
% This function evaluates the conductivity assuming lineare response and a
% broadening parameter Gamma.

N       = 101;                          % Number of eigenstates


% Load Eigenvalues and Dipole Operator
d2      = abs(lho.DipoleOperator).^2;   % Dipole Squared [aL^2]
eng     = lho.EigenValues;              % Energies [Hz]
eng     = eng-eng(1);                   % Offerset energies
eng     = eng(1:N);                     % Energy subspace
d2      = d2(1:N,1:N);                  % Dipole subspace

% Get only specific diagonals (optional)
W = ones(N,N);      % Full weights
S0 = full(spdiags(zeros(N,1),0,W));       % diagonals go to zero; 
S1 = full(spdiags(zeros(N,2),[-1 1],W));  % 1st diagonal to to zero
S2 = full(spdiags(zeros(N,2),[-2 2],W));  % 2nd diagonal to to zero
S3 = full(spdiags(zeros(N,2),[-3 3],W));  % 3rd diagonal to to zero
S4 = full(spdiags(zeros(N,2),[-4 4],W));  % 4th diagonal to to zero
S5 = full(spdiags(zeros(N,2),[-5 5],W));  % 5th diagonal to to zero

% W = logical(~S1);               % Only first 1st diagonal
% W = logical(~S3);               % Only first 3rd diagonal
% W = logical(~S5);               % Only first 5th diagonal
% W = logical(~S3+~S2+~S0);       % Only 0th, 2nd, and 3rd diagonal
% % W = S1;                         % Only non 1st diagonal

% Meshgrid data
[EE1,EE2] = meshgrid(eng,eng);
dEE = EE1-EE2;

% Check for sizing
z = sigma_func(T,G,f);

% Total conductivity wrapper
% We use a wrapper function because f could be a vector
    function y = sigma_func(T,G,f)      
        y = arrayfun(@(f) sigma_helper(T,G,f),f);
    end

% Total conductivity function
    function y = sigma_helper(T,G,f)
        Z = sum(exp(-eng/T),'all');
        A = -1i*f*((exp(-EE1/T)-exp(-EE2/T))/Z).*d2.*W...
            ./((f-dEE)+1i*G/2/(2*pi));
        y = sum(A,'all');
    end

end

