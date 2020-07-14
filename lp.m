function y=lp(nu,x)
    switch nu
        case 0
            y = x.^0;
        case 1
            y = x;
        case 2
            y = 0.5.*(-1 + 3.*x.^2);
        case 3
            y = 0.5.*(-3.*x + 5.*x.^3);
    end
end