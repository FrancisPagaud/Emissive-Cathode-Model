function rep = Ie(n, Te, Vp, Vb, dx)
    Ii = Iis(n,Te,dx);
    len = length(Vb);
    rep = zeros(1,len);
    Mi = 40*1.67e-27;
    me = 9.1e-31;
    Lambda = log(sqrt(2*Mi/(pi*me)));
    
    for i = (1:len)
        if Vb(i) < Vp
            rep(i) = Ii*exp(Lambda + (Vb(i)-Vp)/Te);
        else
            rep(i) = Ii*exp(Lambda);
        end
    end
end