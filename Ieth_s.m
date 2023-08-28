function rep = Ieth_s(Vp, Vb, Te, n, dx)
    len = length(Vb);
    rep = zeros(1,len);
    for i = (1:len)
        if Vb(i) >= Vp 
            rep(i) = 0;
        else
            m_e = 9.1e-31;
            e = 1.602e-19;
            Ce = sqrt(8*e*Te/(pi*m_e));
            r = 2.53e-4;

            phi = (Vb(i)-Vp)/Te;
            F = exp(phi)-1;
            beta_0 = -4*phi.^2-2*phi.*(F.^2-2*F);
            beta_1 = 4*(-2*F-1).*phi.^2 + 8*F.*phi - F.^2;
            beta_2 = 4*phi.^2 - 8*phi.^3;
            G = (-beta_1 + sqrt(beta_1.^2 - 4*beta_0.*beta_2)) ./ (2*beta_2);

            rep(i) = (2*pi*r*dx) * G*e*n*Ce.*sqrt(-pi*phi) ./ (2*(1+G));
        end
    end
end