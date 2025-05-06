function output=bootstrap_linear(x,var,nBootstraps)

x=x(:);
var=var(:);
data=[x var];

[x,inds]=sort(x,'ascend');
var = var(inds);

x0 = median(x);
var0_guess = var(1);
m_guess = (var(end)-var(1))/(x(end)-x(1));


lin_func = @(P,x) ...
    P(1)*(x-x0)+P(2);


P_guess = [m_guess var0_guess];
    
%% Normal Fitting
options = optimset('Display','off');    


[fout,resnorm,residual,exitflag,output0,lambda,jacobian] = ...
    lsqcurvefit(lin_func, P_guess, x,var, [], [], options);
conf = nlparci(fout,residual,'jacobian',jacobian);


m_val = fout(1);
m_err = (conf(1,2)-conf(1,1))/2;

b_val = fout(2);
b_err = (conf(2,2)-conf(2,1))/2;


P_fit=[m_val b_val];
P_err=[m_err b_err];

%% Bootstrap Fitting
    function fittedParams=fitModel(data)
        options = optimset('Display','off');    
        fittedParams = lsqcurvefit(lin_func, P_guess, data(:,1), data(:,2), [], [], options);
    end

% nBootstraps = 1e2;
% Apply bootstrap
[bootstat, bootsam] = bootstrp(nBootstraps, @fitModel, data);


%% Create Ouputs

output = struct;
output.x = x;
output.var = var;
output.FitParam = P_fit;
output.FitErr = P_err;
output.BootStat = bootstat;
output.BootSam = bootsam;
output.FitFunc = lin_func;
end
