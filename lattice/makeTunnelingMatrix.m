function T = makeTunnelingMatrix(n,jmax,tList)
% Kinetic energy tunneling operator with periodic boundary conditions;
%

if length(tList)<jmax
    error('requesting more tunneling elements than you have provided')
end

a = n-1;                        % Useful for later
tList = tList(1:jmax);          % Value up to desired tunneling order
tList = -tList;                 % Kinetic Energy matrix is negative of tunneling
dlist = [flip(tList) tList];    % off diagonal values
L = length(tList);              % Number of tunnelings

diag_list = [-L:-1 1:L];        % Off diagonal list
diag_list_pbc = ...             % Off diagonal list for PBC
    [flip([-a:(-a+(L-1))]) flip((a-(L-1)):a)];

dvecs = repmat(dlist,[n 1]);    % Vectors of all tunneling elements

% Create the Tunneling matrix
T = full(spdiags(dvecs,diag_list,n,n)+spdiags(dvecs,diag_list_pbc,n,n));
end

