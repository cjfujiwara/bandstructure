function l=harmonic_length(U)

out = constants;
m = out.m;
hbar = out.hbar;

omega = 2*pi*sqrt(4*U)*out.fr;

l = 1./sqrt(m*omega/hbar);
end

