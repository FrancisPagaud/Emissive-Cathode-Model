function rep = Ieth(T,dx)
    Ag = 6e5;
    e = 1.602e-19;
    k_b = 1.381e-23;
    W = 4.54;
    r = 2.53e-4;
    rep = Ag*2*pi*r*dx*T.^2.*exp(-(e*W)./(k_b*T));
end