
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
%5
% Er:= hbar^2*k_L^2/(2*m)
%
% In this form the lattice spacing is pi and the FBZ spans [-1,1]

%% Initialize Settings
disp([repmat('-',1,20) ' SETTINGS ' repmat('-',1,20)])
npt=constants;
npt.depth=[2.5 5 7.5 10 15 20]; 
%% Options

opts = struct;
opts.bands = [1 2 3];

%% Calculate Band Structure
disp([repmat('-',1,20) ' BAND STRUCTURE ' repmat('-',1,20)])
tic
close all
opts.doPlot = 0;
[npt,hF1]=getBandStructure(npt,opts);


%% Calculate Tunneling Matrix Elements\
disp([repmat('-',1,20) ' TUNNELING ' repmat('-',1,20)])
% [npt,hF1_a]=calculateAMCoupling(npt,opts);
npt = calculateTunneling(npt);
showTunneling(npt,opts);

%% Interaction
disp([repmat('-',1,20) ' BUSCH ' repmat('-',1,20)])
showTrappedInteraction;

%% Gaps
% [npt,hF2]=calculateGapTunneling(npt);
%% Wannier Calculation
disp([repmat('-',1,20) ' WANNIER ' repmat('-',1,20)])
% Calculate wannier functions from Bloch waves.

opts.doPlot=0;
npt=wannier(npt,opts);

if opts.doPlot;showWannier(npt,opts);end
 
 %% Interaction Integrals
 % Calculate contact interaction hamiltonian
 npt = wannier_overlap_s(npt,opts);
 

 hFa = showSwave(npt,opts);



%% Feshbach Resonance
% Plot the feshbach resonance in free space
% hF4=showSwaveFesh;
% hF7=showPwaveFesh;

% Shift due to zero point energy
% hF = show_lattice_shift_s;

%% Saturation


%% Data RF

% differential_swave_shift2;
% differential_swave_shift3;
% 
% differential_swave_shift4;
% 
% differential_swave_shift5;
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
