function [res] = Integrate(pars); %fids, t1, t2)
    global NMRstat;
    global NMRpar2d;
    global TDfids;
    ints = [];
    
    setup;
    parmaecho (pars{:});
    t2 = NMRstat{2}(1).*1e6;
    t1 = NMRstat{8}.*1e6;
    fids = TDfids;
    
    for i = 1:size(fids,2)
        fid = fids(:,i);
        t = (1:size(fid,1)).*NMRstat{4}.*1e6;
        
        delta_ph = RotatePhase([t' fid],t1,t2);
        ph0 = angle(fid);
        new_fid = abs(fid).*(cos(ph0+delta_ph)+1i*sin(ph0+delta_ph));
        sub_curve = new_fid(find(t >= t1 & t <= t2));
        ints = [ints; trapz(real(sub_curve))];
        figure(2);
        plot(real(new_fid),'Color','red');
        hold on;
        plot(imag(new_fid),'Color','green');
        hold off;
        i
    end
    filename = pars{1};
    prev_data = load(strrep(filename,'.pna','.dat'));
    prev_data = sortrows(prev_data,1);
    tmp = [NMRpar2d.val' ints];
    tmp = sortrows(tmp,1);
    res = [tmp prev_data(:,3)];
end

function [delta_ph] = RotatePhase(fid, t1, t2)
    sub_curve = fid(find(fid(:,1) >= t1 & fid(:,1) <= t2),2);
    ph0 = angle(sub_curve);
    intR_max = trapz(real(sub_curve));
    delta_ph = 0;
    for ph=0:(pi/1000):(2*pi)
        new_curve = abs(sub_curve).*(cos(ph+ph0)+1i*sin(ph+ph0));
        intR = trapz(real(new_curve));
        if intR > intR_max
            intR_max = intR;
            delta_ph = ph;
        end
    end
end