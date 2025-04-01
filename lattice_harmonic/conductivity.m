function conductivity(x,z,src)

% Make into column vector
x=x(:);
z=z(:);

d2 = abs(src.DipoleOperator).^2;
eng = src.EigenValues;
eng = eng-eng(1);

[EE1,EE2] = meshgrid(eng,eng);
dEE = EE1-EE2;

    % function f_boltzman
sR = @(T,G,f) sum(real(-1i*f*(exp(-EE1/T)-exp(-EE2/T)).*d2./((f-dEE)+1i*G/2/(2*pi))),'all');
sI = @(T,G,f) sum(imag(-1i*f*(exp(-EE1/T)-exp(-EE2/T)).*d2./((f-dEE)+1i*G/2/(2*pi))),'all');

    function y = sigma_real(T,G,f)
        y = arrayfun(@(f) sR(T,G,f),f);
    end

    function y = sigma_imag(T,G,f)
        y = arrayfun(@(f) sI(T,G,f),f);
    end

foo = @(P,f) [sigma_real(P(1),P(2),f); sigma_imag(P(1),P(2),f)];

% foo = @(P,f) [arrayfun(@(P,f]
y = [real(z); imag(z)];


P = [500 200];

keyboard
[fout,resnorm,residual,exitflag,output]=lsqcurvefit(foo,P,x,y);

keyboard


% foo=@(P,X) [P(1)*X; P(1)*X.^2];
% Y = [Y1; Y2];

% % fun=@(x,xdata) [f1(xdata(1:numel(xdata)/2),x(1),x(2),x(3),x(4)); f2(xdata(1:numel(xdata)/2),x(1),x(2),x(3),x(4))];
% % x0=[- - - -];
% 





end

