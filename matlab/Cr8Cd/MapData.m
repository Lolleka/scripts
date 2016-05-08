%Crea una mappa 2d dei dati sperimentali, estendendo il risultato
% sperimentale tramite un blur gaussiano 2D dei dati
function [map] = MapData(filename, h, f)
    global bdata;
    %mantiene tutti i dati in una struttura
    %f = frequenza
    %h = campo
    %a = intensità
    bdata = struct('f', [], 'h', [], 'a', [], 'bp', []);
    temp = load(filename)';
    bdata.h = [bdata.h temp(1,:)];
    bdata.f = [bdata.f temp(2,:)];
    bdata.a = [bdata.a temp(3,:)];
    bdata.bp = [bdata.bp temp(4,:)];
    
    %rimuvoe gli elementi con a = 0
    ind = find(bdata.bp == 0);
    
    bdata.a(ind) = [];
    bdata.f(ind) = [];
    bdata.h(ind) = [];
    
    %ordina a in maniera crescente in modo da usare
    % i 4 elementi più piccoli per correggere la baseline
    bdata.a = (bdata.a-min(bdata.a))*((max(bdata.a)-min(bdata.a))^0.2);
    sorted = sort(bdata.a);
    bdata.a = bdata.a - mean(sorted(1:4));
    %azzera gli elementi negativi
    bdata.a = bdata.a .* (bdata.a > 0); 
    
    %inizializza la mappa a 0
    map = zeros(size(h,2),size(f,2));
    
    %parametri del 2D gaussian blur
    sig_f = 1;
    sig_h = 0.1;
    
    %waitbar
    [fp fn fe] = fileparts(filename);
    wb = waitbar (0, ['Blurring data map of ' fn '.' fe ' ...']);
    
    %cicla sulla mappa per riempirla
    for i = 1:size(h,2)
        wb = waitbar (i/size(h,2));
        for j = 1:size(f,2)
           map(i,j) = 0;
           for k = 1:size(bdata.h,2)
               %computa elementi finiti solo entro 4 sigma dalla posizione
               %del dato sperimentale (risparmia tempo)
               if ((abs(h(i)-bdata.h(k)) < 4*sig_h) & (abs(f(j)-bdata.f(k)) < 4*sig_f)) 
                   %asymmetric 2D gaussian blur
                    map(i,j) = map(i,j) + Gauss2D (h(i),f(j),bdata.h(k),bdata.f(k),sig_h,sig_f)*bdata.a(k);
               end
           end
        end
    end
   close(wb);
end

