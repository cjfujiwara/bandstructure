function [wf] = makeBlochState(v,k)
%MAKEWAVEFUNCTION Summary of this function goes here
%   Detailed explanation goes here
    wf=@(x) v(1)*sqrt(1/pi);
    for nn=1:2:length(v)-1
       wf=@(x) wf(x) + ...
           v(nn+1)*sqrt(2/pi)*sin((nn+1).*x)+...
           v(nn+2)*sqrt(2/pi)*cos((nn+1).*x);    
    end
%     wf=@(x) wf(x).*exp(-1i*k.*x*0);

end
