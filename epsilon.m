function rep = epsilon(T)
%     a = 0.105e-3;
%     b = 0.065;
%     rep = a*T+b;
 
%     sigma = 5.67e-8;
%     a = -0.55 + (T-1800)/1000*0.125;        %Pas nécessaire
%     b = -2.91 + (T-1800)/1000*0.73;         %Pas nécessaire
%     lambda = (1e-7:2e-7:1e-4);
%     h = 6.626e-34;
%     c = 2.99792e8;
%     k_b = 1.38e-23;
%     K = 2*pi*h*c^2;
%     e_lambda = a.*log10(lambda)+b;
%     rep = trapz(lambda, K./(lambda.^5).*e_lambda./(exp(h*c./(lambda*k_b*T))-1))/(sigma*T^4);
    
    b1 = -0.1342;       %Fit from Cezairliyan
    b2 = 2.573e-4;
    b3 = -3.481e-8;
    rep = b1 + b2*T + b3*T.^2;
end