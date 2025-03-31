function E = mu2E(E,T,mu)
    E = sum(E./(exp((E-mu)/T)+1));
end

