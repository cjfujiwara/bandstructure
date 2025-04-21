function bs_test
x = (1:1000)';
y1 = 2*x+1 + random('normal',0,2,1000,1);
y2 = 0.5*x - 2 + random('normal',0,2,1000,1);
z = 3 * y1 + 2 * y2 + 1 + random('normal',0,2,1000,1);
% look at the very different outputs of regress vs fitlm:
b11 = regress(z,[y1 y2 ones(1000,1)])
b21 = fitlm([y1 y2],z)


% One way around the problem is to define & bootstrap your own version
% of the function fitlm that returns only a vector of values
b2 = bootstrp(1000, @myfitlm, [y1 y2], z); % < -- This works
function myout = myfitlm(preds,tobepred)
% select out the coefficients that are to be bootstrapped.
b21 = fitlm(preds,tobepred);
myout = b21.Coefficients.Estimate;
end
end

