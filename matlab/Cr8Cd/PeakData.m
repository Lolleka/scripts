%trova i picchi nei file di dati
function [mat_res] = PeakData(filename, f_or_h)
    global bdata;
    %mantiene tutti i dati in una struttura
    %f = frequenza
    %h = campo
    %a = intensità
    bdata = struct('f', [], 'h', [], 'a', []);
    temp = load(filename)';
    bdata.h = [bdata.h temp(1,:)];
    bdata.f = [bdata.f temp(2,:)];
    bdata.a = [bdata.a temp(3,:)];
    
    %rimuvoe gli elementi con a = 0
    ind = find(bdata.a == 0);
    bdata.a(ind) = [];
    bdata.f(ind) = [];
    bdata.h(ind) = [];
    
    %ordina a in maniera crescente in modo da usare
    % i 4 elementi più piccoli per correggere la baseline
    sorted = sort(bdata.a);
    bdata.a = bdata.a - mean(sorted(1:4));
    %azzera gli elementi negativi
    bdata.a = bdata.a .* (bdata.a > 0); 
    
    N_interp = 100;
    if f_or_h == 'f'
        min_f = min(bdata.f);
        max_f = max(bdata.f);
        ndata = struct('f', [], 'h', [], 'a', []);
        ndata.f = min_f:((max_f-min_f)/(N_interp-1)):max_f;
        ndata.a = interp1(bdata.f,bdata.a,ndata.f)./max(bdata.a);
        ndata.h = interp1(bdata.f,bdata.h,ndata.f);
        %figure(1);
        %plot(ndata.f,ndata.a,'o');
        [pks,pos_f]=findpeaks(ndata.a,'MINPEAKHEIGHT',0.1,'MINPEAKDISTANCE',6);
        pos_h = repmat(ndata.h(1),1,size(pos_f,2));
        pos_f = ndata.f(pos_f);
        %[pks,locs]=findpeaks(ndata.a);
        %hold on;
        %p = plot(ndata.f(locs),pks,'o');
        %set (p,'Color','red');
        %hold off;
    elseif f_or_h == 'h'
        min_h = min(bdata.h);
        max_h = max(bdata.h);
        ndata = struct('f', [], 'h', [], 'a', []);
        ndata.h = min_h:((max_h-min_h)/(N_interp-1)):max_h;
        ndata.a = interp1(bdata.h,bdata.a,ndata.h)./max(bdata.a);
        ndata.f = interp1(bdata.h,bdata.f,ndata.h);
        %figure(1);
        %plot(ndata.f,ndata.a,'o');
        [pks,pos_h]=findpeaks(ndata.a,'MINPEAKHEIGHT',0.1,'MINPEAKDISTANCE',6);
        pos_f = repmat(ndata.f(1),1,size(pos_h,2));
        pos_h = ndata.h(pos_h);
        %[pks,locs]=findpeaks(ndata.a);
        %hold on;
        %p = plot(ndata.h(pos_h),pks,'o');
        %set (p,'Color','red');
        %hold off;        
    end
    mat_res = [pos_h' pos_f' pks'];
end