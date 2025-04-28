function bootstrap_com(digdata)

    function out=mean2D(data)
        % options = optimset('Display','off');    
        % fittedParams = lsqcurvefit(oscillations_wrapper, P_guess, data(:,1), data(:,2), [], [], options);
        out = [mean(data(:,1)) mean(data(:,2))];
    end
    freqs=[];
    V2=[];
for nn=1:length(digdata)
    Xc=[];
    Yc=[];
    for ii=1:length(digdata(nn).Ratom)
        Rimg=digdata(nn).Ratom{ii};
        X = Rimg(1,:);X=X(:);
        Y = Rimg(2,:);Y=Y(:);
        data=[X Y];
                
        [bootstat,bootsam] = bootstrp(100,@mean2D,data);
        pdXc = fitdist(bootstat(:,1),'Normal');
        pdYc = fitdist(bootstat(:,2),'Normal');
        px = paramci(pdXc);
        py = paramci(pdYc);

        Xc(ii) = pdXc.mu;
        Yc(ii) = pdYc.mu;
        XcErr(ii) = (px(2,1)-px(1,1))*0.5;
        YcErr(ii) = (py(2,1)-py(1,1))*0.5;
    end
    Xc=Xc(:);
    XcErr=XcErr(:);
    Yc=Yc(:);
    YcErr=YcErr(:);

    P = [digdata(nn).Params];
    freq = P(1).conductivity_mod_freq;

    px_per_site = 2.68;
    um_per_site = 0.532;
    um_per_px   = um_per_site/px_per_site;

    T = [P.conductivity_mod_time];
    Tr = [P.conductivity_mod_ramp_time];
    Ttot = T+Tr;
    Ttot=Ttot(:);
    freqs(nn) = freq;
    output(nn)=bootstrap_oscillations_fit(Ttot,Xc*um_per_px,Yc*um_per_px,freq*1e-3);
    V2(nn) = [P(1).conductivity_ODT2_mod_amp];
end
allP = [output.FitParam];

S = allP(1:4:end);
C = allP(2:4:end);
keyboard

end

