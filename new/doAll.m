
disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))    

%% Comment on Units
% This code calculates values paramtrized by the lattice momentum hbar*k_L
% so that adopt the dimensionless position and momentum 
%
% X:=k_L*x, P:=p/(hbar*k_L) ==> [X,P]=i
%
% So that the dimensionless Hamiltonian takes the form
%
% h = U_0\cos^2(X)+P^2
%
% Where the lattice depth U_0 is parametrized with the recoil energy
%
% Er:= hbar^2*k_L^2/(2*m)
%
% In this form the lattice spacing is pi and the FBZ spans [-1,1]
%%

close all


npt=constants;
npt.depth=200;

[npt,hF1]=getBandStructure(npt);
[npt,hF2]=calculateGapTunneling(npt);
[npt,hF3]=wannier(npt);

%% Overlap Integrals

npt = wannier_overlap(npt);

%% Feshbach Resonance
% Plot the feshbach resonance in free space
hF4=showSwave;
hF7=showPwave;

% Shift due to zero point energy
hF = show_lattice_shift_s;

%% Saturation

hF5=showTrappedInteraction;

%% Data RF

differential_swave_shift;

%% Data Ideal
s_wave_shift;

%%
doSave=0;
if doSave
    disp('saving figures ...');    
    figs=get(groot,'children');
    
    for kk=1:length(figs)
        fprintf(['saving figure ' num2str(kk) ' ... ']);
        print(figs(kk),['figs/' figs(kk).Name],'-dpng','-r400'); 
        savefig(figs(kk),['figs/' figs(kk).Name]); 
        disp('done');

    end    
end
