function [res] = EprGetH(nu,g)
    h = 6.62e-34;
    mub = 9.27E-24;
    S = 1/2;
    res = h*nu / (2*S*mub*g);
end