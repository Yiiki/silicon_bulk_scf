function y=bj(nu,x)
    switch nu
        case 0
            y = sin(abs(x))./abs(x);
        case 1
            y = 1./x.*(-cos(abs(x)) + sin(abs(x))./abs(x));
        case 2
            y = 1./abs(x) .* (-3.*cos(abs(x))./abs(x) - sin(abs(x))+3.*sin(abs(x))./(abs(x)).^2);
    end
end