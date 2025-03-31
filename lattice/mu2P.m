function P = mu2P(E,T,mu)
    P = 1./(exp((E-mu)/T)+1);
end