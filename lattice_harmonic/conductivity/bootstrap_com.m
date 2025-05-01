function out=bootstrap_com(digdata)
% Image Calibrations
px_per_site = 2.68;
um_per_site = 0.532;
um_per_px   = um_per_site/px_per_site;


disp('Bootstrap Conductivity')
freqs = zeros(length(digdata),1);
V2=[];

N_boot_moments      = 1000;
N_boot_oscillations = 1000;

CoM_range           = 4;

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

out = struct;

% Iterate over all digdata
for nn=1:length(digdata)
    disp([num2str(nn) ' of ' num2str(length(digdata))])    

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

     for ii=1:length(digdata(nn).Ratom)
         fprintf(['image ' num2str(ii) '/' num2str(length(digdata(nn).Ratom))]);
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
        disp(' done');
     end
    fprintf('bootstrapping heating ...');
    heatX(nn)=bootstrap_linear(Ttot,Xvar);
    heatY(nn)=bootstrap_linear(Ttot,Yvar);
    disp('done');


    % heatXY(nn) = bootstrap_linear2d(Ttot,[Xvar Yvar]); 
    fprintf('bootstrapping oscillations ...');
    output(nn)=bootstrap_oscillations_fit(Ttot,Xcom,freq*1e-3);
    disp('done');
    C=output(nn).Covariance; % covariance matrix
    

    
    % Atom Number and Gauss Charge Density
    N =[digdata(nn).Natoms];
    sx = sqrt(Xvar);
    sy = sqrt(Yvar);
    
    rho_charge_gauss_peak = N./(2*pi*sx.*sy);
%     rho_charge_gauss_avg  = N./
    
    % 2D gaussform N/(sqrt(2*pi*sx^2)*sqrt(2*pi*sx^2)) exp(-x^2/2s^2)exp(-y^2/2s^2)

    

    hF = figure(nn);
    clf
    set(hF,'color','w');

    tg = uitabgroup(hF);
    summary=uitab(tg,'Title','summary','backgroundcolor','w');
    x_osc=uitab(tg,'Title','oscillations','backgroundcolor','w');
    t_heating=uitab(tg,'Title','heating','backgroundcolor','w');


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




    tt=linspace(min(Ttot),max(Ttot),100);
    plot(tt,output(nn).FitFunc(output(nn).FitParam,tt),'k-','parent',ax_comx,'linewidth',1);

    plot(tt,heatX(nn).FitFunc(heatX(nn).FitParam,tt),'k-','parent',ax_varx,'linewidth',1);
    plot(tt,heatY(nn).FitFunc(heatY(nn).FitParam,tt),'k-','parent',ax_vary,'linewidth',1);


    subplot(2,4,1,'parent',x_osc)
    histogram(output(nn).BootStat(:,1));
    xlabel('S (um)')
    title('com-displacement in-phase ')

    subplot(2,4,2,'parent',x_osc)
    histogram(output(nn).BootStat(:,2));
    xlabel('C (um)')
    title('com-displacement out-phase ')

    subplot(2,4,3,'parent',x_osc)
    histogram(output(nn).BootStat(:,3));
    xlabel('center (um)')
    title('center position')

    subplot(2,4,4,'parent',x_osc)
    histogram(output(nn).BootStat(:,4));
    xlabel('velocity (um/ms)')
    title('u-scope drift "velocity"')

    subplot(2,4,5,'parent',x_osc)
    plot(output(nn).BootStat(:,1),output(nn).BootStat(:,2),'.');
    xlabel('S (\mum)')
    ylabel('C (\mum)')
    title('S-C covariance')
    axis equal tight
    hold on

    c1=mean(output(nn).BootStat(:,1));
    c2=mean(output(nn).BootStat(:,2));

    xa = c1 + [-1 1]*C(1,1);
    ya = c2 +[-1 1]*C(1,2);

    xb = c1 + [-1 1]*C(1,2);
    yb = c2 +[-1 1]*C(2,2);


    subplot(2,4,6,'parent',x_osc)
    plot(output(nn).BootStat(:,3),output(nn).BootStat(:,4),'.');
    xlabel('center')
    ylabel('velocity')
    title('x0-v covariance')


    subplot(2,4,[7 8],'parent',x_osc)
    sC=num2str(round(C,4));
    text(.5,.5,sC,'units','normalized',...
        'horizontalalignment','center','verticalalignment','middle')
    title('S-C-x0-v0 covariance matrix')

    subplot(2,3,1,'parent',t_heating)
    histogram(heatX(nn).BootStat(:,1))
    xlabel('slope um^2/ms')
    title('X var heating')

    subplot(2,3,4,'parent',t_heating)
    histogram(heatX(nn).BootStat(:,2))
    xlabel('offset um^2')
    title('X var offset')


    subplot(2,3,2,'parent',t_heating)
    histogram(heatY(nn).BootStat(:,1))  
    xlabel('slope um^2/ms')
    title('Y var heating')


    subplot(2,3,5,'parent',t_heating)
    histogram(heatY(nn).BootStat(:,2))  
    histogram(heatX(nn).BootStat(:,2))
    xlabel('offset um^2')
    title('Y var offset')

    subplot(2,3,3,'parent',t_heating)
    % plot(heatX(nn).BootStat(:,1),heatY(nn).BootStat(:,1),'o')
    xlabel('x slope um^2')
    ylabel('y slope um^2')
    axis equal tight
    title('X-Y heating correlator')
    
    out(nn).SourceDirectory = digdata(nn).SourceDirectory;
    out(nn).FileNames = digdata(nn).FileNames;
    out(nn).Params = digdata(nn).Params;
    out(nn).Freqency  = freq;
    out(nn).S = output(nn).FitParam(1);
    out(nn).SErr = output(nn).FitErr(1);
    out(nn).C = output(nn).FitParam(2);
    out(nn).CErr = output(nn).FitErr(2);
    out(nn).BootStrapOscillations = output(nn);
    out(nn).BootStrapHeatX = heatX(nn);
    out(nn).BootStrapHeatY = heatY(nn);
% nn
% keyboard
end
allP = [output.FitParam];
 
S = allP(1:4:end);
C = allP(2:4:end);
% keyboard

end

