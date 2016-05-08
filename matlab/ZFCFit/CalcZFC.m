function [ res ] = CalcZFC(T, A, nu0, nuT, Kg, KDelta)
    kB = 1.3806503E-23;
    K = 0.1:0.1:10000;    
    %distr = (1 ./ (K .* sigma * sqrt(2 * pi))) .* exp(-((log(K) - log(Kg)) .^ 2) ./ (2 * sigma * sigma));
    distr = my_distr(K,Kg,KDelta);
    
    ZFCSPM = K.^2 .* distr;
    ZFCBlocked = K .* distr;
    KLim = (0.9609 .* log(nu0 .* T ./ nuT) - 1.629) .* T;
    
    res_ZFCSPM = zeros(1,length(T));
    res_ZFCBlocked = zeros(1,length(T));
    res_ZFCSPM_tot = trapz(K, ZFCSPM);
    for i = 1:length(T)
        ind_a = find(K<=KLim(i));
        ind_b = find(K>KLim(i));
        
        if not(isnan(ind_a))
            if KLim(i) <= K(length(K))
                newK_a = [K(ind_a) KLim(i)];
                newZFCSPM_point = interp1(K, ZFCSPM, KLim(i));
                newZFCSPM = [ZFCSPM(ind_a) newZFCSPM_point];
                res_ZFCSPM(i) = trapz(newK_a, newZFCSPM);
            else
                res_ZFCSPM(i) = res_ZFCSPM_tot;
            end
        else
            res_ZFCSPM(i) = 0;
        end
       
        if not(isnan(ind_b))
            newK_b = [KLim(i) K(ind_b)];
            newZFCBlocked_point = interp1(K, ZFCBlocked, KLim(i));
            newZFCBlocked = [newZFCBlocked_point ZFCBlocked(ind_b)];
            res_ZFCBlocked(i) = trapz(newK_b, newZFCBlocked);
        else
            res_ZFCBlocked(i) = 0;
        end        
    end
    res = A * (res_ZFCSPM ./ (kB .* T) + res_ZFCBlocked);
end

