
disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))    

%%

close all

npt=constants;
npt.depth=200;

[npt,hF1]=getBandStructure(npt);
[npt,hF2]=calculateGapTunneling(npt);
[npt,hF3]=wannier(npt);

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
