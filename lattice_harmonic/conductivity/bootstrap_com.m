function [output,hF]=bootstrap_com(composite_data)

digdata=[composite_data.digdata];
% Image Calibrations
px_per_site = 2.68;
um_per_site = 0.532;
um_per_px   = um_per_site/px_per_site;

disp('Bootstrap Conductivity')
freqs = zeros(length(digdata),1);
V2=[];

N_boot_moments      = 100;
N_boot_heating      = 100;
N_boot_oscillations = 100;
CoM_range           = 4;


%% Moment Functions
% centre-of-mass
    function out = mean2D(data)
        out = [mean(data(:,1)) mean(data(:,2))];
    end

% variance
    function out = var2D(data,mu)
        x2 = (data(:,1)-mu(1)).^2;
        y2 = (data(:,2)-mu(2)).^2;
        r2 = x2+y2;
        out = [mean(x2) mean(y2) mean(r2)];
    end

% skew
    function out = skew2D(data,mu,sigma)
        x3 = ((data(:,1)-mu(1))./sigma(1)).^3;
        y3 = ((data(:,2)-mu(2))./sigma(2)).^3;
        out = [mean(x3) mean(y3)];
    end
%% Init Graphics

%%
output = struct;

% Iterate over all composite_data
for jj=1:length(composite_data)
    % Initialize figure for specific spectrum
    disp(['composite_data ' num2str(jj) '/' num2str(length(composite_data))]);
    hF(jj) = figure(200+jj);
    hF(jj).Color='w';
    hF(jj).Name = [composite_data(jj).Name ' moments'];
    hF(jj).WindowStyle = 'docked';
    clf(hF(jj))   
    tab_group_spectrum = uitabgroup(hF(jj));
    drawnow;
    % Get specific digdata
    digdata=[composite_data(jj).digdata];
    for nn=1:length(digdata)        
        freq_drive = [digdata(nn).Params(1).conductivity_mod_freq];
        output(jj).Frequency_Hz(nn) = freq_drive;

        fprintf([num2str(freq_drive) ' Hz:']);

        % Initialize Panel Objects
        panel_frequency = uitab(tab_group_spectrum,'title',[num2str(freq_drive) ' Hz']);
        tab_group_frequency = uitabgroup(panel_frequency);
        summary=uitab(tab_group_frequency,'Title','summary','backgroundcolor','w');
        tab_osc=uitab(tab_group_frequency,'Title','oscillations','backgroundcolor','w');
        tab_heating=uitab(tab_group_frequency,'Title','heating','backgroundcolor','w');    
        
        % First Moment Initialize
        Xcom      = zeros(size(digdata(nn).Ratom,2),1);
        XcomErr   = zeros(size(digdata(nn).Ratom,2),1);
        Ycom      = zeros(size(digdata(nn).Ratom,2),1);
        YcomErr   = zeros(size(digdata(nn).Ratom,2),1);
    
        % Second Moment Initialize
        Xvar      = zeros(size(digdata(nn).Ratom,2),1);
        XvarErr   = zeros(size(digdata(nn).Ratom,2),1);
        Yvar      = zeros(size(digdata(nn).Ratom,2),1);
        YvarErr   = zeros(size(digdata(nn).Ratom,2),1);
        Rvar      = zeros(size(digdata(nn).Ratom,2),1);
        RvarErr   = zeros(size(digdata(nn).Ratom,2),1);
    
        % Third Moment Initialize
        Xskew     = zeros(size(digdata(nn).Ratom,2),1);
        XskewErr  = zeros(size(digdata(nn).Ratom,2),1);
        Yskew     = zeros(size(digdata(nn).Ratom,2),1);
        YskewErr  = zeros(size(digdata(nn).Ratom,2),1);
    
        % Get relevant data
        P = [digdata(nn).Params];
        freq = P(1).conductivity_mod_freq;
        T       = [P.conductivity_mod_time];
        Tr      = [P.conductivity_mod_ramp_time];
        Ttot    = T+Tr;
        Ttot    = Ttot(:);
        freqs(nn) = freq;
    
          
        fprintf('bootstraping... ')
         fprintf('moments... ')

         for ii=1:length(digdata(nn).Ratom)
            data = digdata(nn).Ratom{ii}'*um_per_px;            % data             
    
             % Bootstrap centre-of-mass
            [bootstat_com,bootsam_com] = ...
                bootstrp(N_boot_moments,@mean2D,data);   % bootstrap
            pdXcom = fitdist(bootstat_com(:,1),'Normal');       % x gauss stats
            pdYcom = fitdist(bootstat_com(:,2),'Normal');       % y gauss stats
            px = paramci(pdXcom);                               % x gauss conf
            py = paramci(pdYcom);                               % y gauss conf
            Xcom(ii) = pdXcom.mu;                               % x com
            Ycom(ii) = pdYcom.mu;                               % y com
            XcomErr(ii) = (px(2,1)-px(1,1))*0.5;                % x com err
            YcomErr(ii) = (py(2,1)-py(1,1))*0.5;                % y com err
    
            % Bootstrap the variance
            [bootstat_var,bootsam_var] = ...
                bootstrp(N_boot_moments,@(x) var2D(x,[Xcom(ii) Ycom(ii)]),data);
            pdXvar = fitdist(bootstat_var(:,1),'Normal');
            pdYvar = fitdist(bootstat_var(:,2),'Normal');
            pdRvar = fitdist(bootstat_var(:,3),'Normal');
            px = paramci(pdXvar);                               % x gauss conf
            py = paramci(pdYvar);                               % y gauss conf
            pr = paramci(pdRvar);                               % y gauss conf
    
            Xvar(ii) = pdXvar.mu;                               % x var
            Yvar(ii) = pdYvar.mu;                               % y var
            Rvar(ii) = pdRvar.mu;                               % r var
            XvarErr(ii) = (px(2,1)-px(1,1))*0.5;                % x var err
            YvarErr(ii) = (py(2,1)-py(1,1))*0.5;                % y var err
            RvarErr(ii) = (pr(2,1)-pr(1,1))*0.5;                % y var err    
    
            % Bootstrap the Skew
            [bootstat_skew,bootsam_skew] = ...
                bootstrp(N_boot_moments,...
                @(x) skew2D(x,[Xcom(ii) Ycom(ii)],[sqrt(Xvar(ii)) sqrt(Yvar(ii))]),...
                data);
            pdXskew = fitdist(bootstat_skew(:,1),'Normal');
            pdYskew = fitdist(bootstat_skew(:,2),'Normal');
            px = paramci(pdXskew);                               % x gauss conf
            py = paramci(pdYskew);                               % y gauss conf
            Xskew(ii) = pdXskew.mu;                              % x skew
            Yskew(ii) = pdYskew.mu;                              % y skew
            XskewErr(ii) = (px(2,1)-px(1,1))*0.5;                % x skew err
            YskewErr(ii) = (py(2,1)-py(1,1))*0.5;                % y skew err
         end

         fprintf('heating... ')
        bs_heatX=bootstrap_linear(Ttot,Xvar,N_boot_heating);
        bs_heatY=bootstrap_linear(Ttot,Yvar,N_boot_heating);
    
         fprintf('oscillations... ')
        bs_osc=bootstrap_oscillations_fit(Ttot,Xcom,freq*1e-3,N_boot_oscillations);
        C=bs_osc.Covariance; % covariance matrix      
        
        % Atom Number and Gauss Charge Density
        N =[digdata(nn).Natoms];N=N(:);
        sx = sqrt(Xvar);sx=sx(:);
        sy = sqrt(Yvar);sy=sy(:);    
        rho_charge_gauss_peak = N./(2*pi*(sx/um_per_site).*(sy/um_per_site));       

        %% Assign Output
        
        output(jj).SourceDirectory{nn} = {digdata(nn).SourceDirectory};
        output(jj).FileNames{nn} = {digdata(nn).FileNames};
        output(jj).Params{nn} = [digdata(nn).Params];
        % output(jj).Depth_Er(nn) = unique([bs_moments(1).Params{1}.lattice_load_depthX]);
        % output(jj).Field_Gauss = unique([bs_moments(1).Params{1}.conductivity_FB_field])+0.123;
        output(jj).N{nn}= N;
        output(jj).Density_PeakGayssCharge{nn} = rho_charge_gauss_peak;
        output(jj).CentreX{nn}                  = Xcom;
        output(jj).CentreY{nn}                  = Ycom;
        output(jj).VarianceX_um{nn}             = Xvar;
        output(jj).VarianceY_um{nn}             = Yvar;
        output(jj).SkewX{nn}                    = Xskew;
        output(jj).SkewY{nn}                    = Yskew;
        output(jj).S_um(nn)                     = bs_osc.FitParam(1);
        output(jj).SErr_um(nn)                  = bs_osc.FitErr(1);
        output(jj).C_um(nn)                     = bs_osc.FitParam(2);
        output(jj).CErr_um(nn)                  = bs_osc.FitErr(2);
        output(jj).HeatX_um2perms(nn)           = bs_heatX.FitParam(1);
        output(jj).HeatXErr_um2perms(nn)        = bs_heatX.FitErr(1);
        output(jj).HeatY_um2perms(nn)           = bs_heatY.FitParam(1);
        output(jj).HeatYErr_um2perms(nn)        = bs_heatY.FitErr(1);

        fprintf(' plotting')
        %% Plot Moments
        % Center of Mass X
        ax_comx=subplot(4,2,1,'parent',summary);
        co=get(gca,'colororder');
        errorbar(Ttot,Xcom,XcomErr,'o','markerfacecolor',co(1,:),...
            'markeredgecolor',co(1,:)*.5,'linewidth',1,'markersize',8,...
            'color',co(1,:)*.5);
        ylim(mean(Xcom)+[-.5 .5]*CoM_range)
        title('X center-of-mass')
        ylabel('E[x] (\mum)')
        hold on
    
        % Center of Mass Y
        subplot(4,2,2,'parent',summary)
        errorbar(Ttot,Ycom,YcomErr,'o','markerfacecolor',co(2,:),...
            'markeredgecolor',co(2,:)*.5,'linewidth',1,'markersize',8,...
            'color',co(2,:)*.5);
        ylim(mean(Ycom)+[-.5 .5]*CoM_range)
        title('Y center-of-mass')
        ylabel('E[y] (\mum)')
    
        % Variance X
        ax_varx=subplot(4,2,3,'parent',summary);
        errorbar(Ttot,Xvar,XvarErr,'o','markerfacecolor',co(1,:),...
            'markeredgecolor',co(1,:)*.5,'linewidth',1,'markersize',8,...
            'color',co(1,:)*.5);
        ylabel('E[(x-\mu)^2] (\mum^2)')
        title('X variance')
        hold on
    
        % Variance Y
        ax_vary=subplot(4,2,4,'parent',summary);
        errorbar(Ttot,Yvar,YvarErr,'o','markerfacecolor',co(2,:),...
            'markeredgecolor',co(2,:)*.5,'linewidth',1,'markersize',8,...
            'color',co(2,:)*.5);
        ylabel('E[(y-\mu)^2] (\mum^2)')
        title('Y variance')
        hold on
    
        % Skew X
        subplot(4,2,5,'parent',summary)
        errorbar(Ttot,Xskew,XskewErr,'o','markerfacecolor',co(1,:),...
            'markeredgecolor',co(1,:)*.5,'linewidth',1,'markersize',8,...
            'color',co(1,:)*.5);
        title('X skew')
        ylabel('E[(x-\mu)^3/\sigma^3]')
    
        % Skew Y
        subplot(4,2,6,'parent',summary)
        errorbar(Ttot,Yskew,YskewErr,'o','markerfacecolor',co(2,:),...
            'markeredgecolor',co(2,:)*.5,'linewidth',1,'markersize',8,...
            'color',co(2,:)*.5);
        title('Y skew')
        ylabel('E[(x-\mu)^3/\sigma^3]')
    
        % Atom Number
        subplot(4,2,7,'parent',summary)
        plot(Ttot,[digdata(nn).Natoms],'o','markerfacecolor',[.5 .5 .5],...
            'markeredgecolor','k','linewidth',1,'markersize',8,...
            'color','k');
        ylabel('number atoms')
        title('atom number');
        ylim([0 max([digdata(nn).Natoms])*1.2]);

        % Density
        subplot(4,2,8,'parent',summary)
        plot(Ttot,rho_charge_gauss_peak,'o','markerfacecolor',[.5 .5 .5],...
            'markeredgecolor','k','linewidth',1,'markersize',8,...
            'color','k');
        ylabel('peak charge density - gauss 2D')
        title('density');
        % ylim([0 max([digdata(nn).Natoms])*1.2]);       
    
        %% Plot Oscillations    
        tt=linspace(min(Ttot),max(Ttot),100);
        plot(tt,bs_osc.FitFunc(bs_osc.FitParam,tt),'k-','parent',ax_comx,'linewidth',1);    
        plot(tt,bs_heatX.FitFunc(bs_heatX.FitParam,tt),'k-','parent',ax_varx,'linewidth',1);
        plot(tt,bs_heatY.FitFunc(bs_heatY.FitParam,tt),'k-','parent',ax_vary,'linewidth',1);        
        subplot(2,4,1,'parent',tab_osc)
        histogram(bs_osc.BootStat(:,1));
        xlabel('S (um)')
        title('com-displacement in-phase ')    
        subplot(2,4,2,'parent',tab_osc)
        histogram(bs_osc.BootStat(:,2));
        xlabel('C (um)')
        title('com-displacement out-phase ')    
        subplot(2,4,3,'parent',tab_osc)
        histogram(bs_osc.BootStat(:,3));
        xlabel('center (um)')
        title('center position')    
        subplot(2,4,4,'parent',tab_osc)
        histogram(bs_osc.BootStat(:,4));
        xlabel('velocity (um/ms)')
        title('u-scope drift "velocity"')    
        subplot(2,4,5,'parent',tab_osc)
        plot(bs_osc.BootStat(:,1),bs_osc.BootStat(:,2),'.');
        xlabel('S (\mum)')
        ylabel('C (\mum)')
        title('S-C covariance')
        axis equal tight
        hold on
        subplot(2,4,6,'parent',tab_osc)
        plot(bs_osc.BootStat(:,3),bs_osc.BootStat(:,4),'.');
        xlabel('center')
        ylabel('velocity')
        title('x0-v covariance')       
        subplot(2,4,[7 8],'parent',tab_osc)
        sC=num2str(round(C,4));
        text(.5,.5,sC,'units','normalized',...
            'horizontalalignment','center','verticalalignment','middle')
        title('S-C-x0-v0 covariance matrix')

        %% Plot Heating   
        subplot(2,3,1,'parent',tab_heating)
        histogram(bs_heatX.BootStat(:,1))
        xlabel('slope um^2/ms')
        title('X var heating')    
        subplot(2,3,4,'parent',tab_heating)
        histogram(bs_heatX.BootStat(:,2))
        xlabel('offset um^2')
        title('X var offset') 
        subplot(2,3,2,'parent',tab_heating)
        histogram(bs_heatY.BootStat(:,1))  
        xlabel('slope um^2/ms')
        title('Y var heating')   
        subplot(2,3,5,'parent',tab_heating)
        histogram(bs_heatY.BootStat(:,2))  
        xlabel('offset um^2')
        title('Y var offset')        
        drawnow;
        disp(' done');
    end
end
% allP = [output.FitParam];
% 
% S = allP(1:4:end);
% C = allP(2:4:end);
% % keyboard

end

