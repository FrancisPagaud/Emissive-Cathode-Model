function rep = Iis(n, Te,dx)
    alpha = 1;
    r = 2.53e-4;
    e = 1.602e-19;
    Mi = 40*1.67e-27;
    u_b = sqrt(e*Te/Mi);
    
    rep = alpha*n*e*(2*pi*r*dx).*u_b;
end