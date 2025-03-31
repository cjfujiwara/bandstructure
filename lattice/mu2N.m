
function Nme = mu2N(E,T,mu)
    Nme = sum(1./(exp((E-mu)/T)+1));
end
