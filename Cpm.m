function rep = Cpm(T)
    rho_m = 19.3e6 ; 
    rep = (0.2294 - 8.24e-5*T + 2.704e-8*T.^2)*rho_m;
end