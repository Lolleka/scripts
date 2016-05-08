function [res] = my_distr(x, dmean, dwidth)
    a = exp(2*log(dmean));
    dsigma = sqrt(log((a + sqrt(a^2+4*a*(dwidth^2))) / (2*a)));    
    res = (1 ./ (x .* dsigma * sqrt(2 * pi))) .* exp(-((log(x) - log(dmean)) .^ 2) ./ (2 * dsigma * dsigma));
end

function [res] = my_distr2(x, dmean, dwidth)
    res = gaussmf(x,[dwidth dmean]);
end