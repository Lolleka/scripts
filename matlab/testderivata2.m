function testderivata (filename)    %evoluzione dello script "script2", i.e. creo una funzione che chiamo "testderivata"
                                    %che vado ad applicare ai file di dati
                                    %con il loro nome vero ossia "filename"                              
    % preparo il file di dati x essere macinato da matlab, i.e.
    % tolgo un invio dalla prima riga e sostituisco la stringa letterale con uno 0 a inizio file
    % rinomino il file cos� modificato "matrix.dat"
    % questo file dati verr� richiamato dallo script
    % perci� posso si modificarlo questo nome, ma allora devo modificarlo tutte le volte che compare nel restro dello script
  
    close all;
    
    punti = find(filename=='.');
    if length(punti)>0
        radice = filename(1:(min(punti)-1));
    else
        radice = filename;
    end
    
    matrix = load(filename);    % battezzo "matrix" il filename che carico
    [row,col]=size(matrix);     % recupera le dimensioni della matrice del file dati ribattezzato "matrix"
                                % e inserisce questi due numeri (righe e colonne) nelle variabili row e col
    x=matrix(1,2:col);
    y=matrix(2:row,1);          % recuperano rispettivamente la prima riga e la prima colonna del mio file dati
                                % tralasciando il primo elemento, che ho settato a 0

    z=matrix(2:row,2:col);      % recupera il resto del mio file dati e lo inserisce in una matrice di nome z
    %figure(1);                  % seleziona una finestra per fare la prima figura
    %imagesc(x,y,z);             % usa i vettori x ed y per stampare in figura la mappa contenuta in z, producendo un plot a colori come in  origin
    %caxis([-0.1 0.1]);          % definisce il range della scala di colore   
    %colorbar
    %xlim([-1 1])                % setta il range dell'asse x
    
    d=[];           % non � inutile, battezza una nuova matrice di nome d, vuota
    row=row-1;
    col=col-1;      % diminuisce di 1 il numero di colonne e righe per le operazioni successive, 
                    % perch� la prima colonna e la prima riga del file dati sono state utilizzate per i vettori x ed y
                    % e non rientrano tra i dati di cui calcolare la derivata

     % il ciclo for successivo elabora tutte le righe della matrice z e fa le seguenti 2 cose:
     % calcola la derivata parziale di z rispetto a x con diff(z(i,1:col))./diff(x)
     %e assegna il risultato a un vettore di nome u, il cui primo elemento � 0
     % poi concatena questo vettore u alla matrice d, in modo da crearla riga per riga
     % cos� ho costruito la matrice contenente le derivate parziali dei miei dati, sulle righe               

    for i=1:row
        u=[0 diff(z(i,1:col))./diff(x)];
        d=[d; u];
    end
    %fighandle = figure(1);          % seleziona un secondo slot per plottare la matrice d
    %N=size(d(i,1:col),2);
    
    %imagesc(x,y,d);     % plotta la matrice d, ovvero la conduttanza differenziale
    %colorbar
    %caxis([0 0.2]);      
    %xlim([-1 1])
    
    %saveas(1,[radice '_plot.jpg']);     % salva la figura 1 in formato jpg, col nome del file cio� "radice" seguito da "_plot"
    %saveas(2,[radice '_der.jpg']);      % salva la figura 2 in formato jpg, col nome del file cio� "radice" seguito da "_der"
    
    %clear d;
    N=size(z(i,1:col),2);
    M=size(z,1);
    k=(1:N)-ceil(N/2);
    k1=0.1*N/2;
    k1=35;
    a=1/(k1^2);
    %par=1-max(0,-a*((abs(x)-40).*(abs(x)-40))+1);
    
    %par=max(0,-a*(x).*(x)+1);
    par=(k>-k1 & k<k1);
    %if abs(k)<k1
    %    par=1;
    %else  
    %    par=0;
    %end
    
    %figure(1);
    %plot(k,par);
    
    f=[];
    fs=[];
    for i=1:row
        unfiltered=fft(z(i,1:N));
        filtered=fftshift(unfiltered).*par;
        f=[f; real(ifft(ifftshift(filtered)))];
        fs=[fs; real(filtered)];
    end
    figure(2);
    
    imagesc(x,y,f);
    %caxis([-0.2 0.2]);
    %colorbar;
    
    fighandle = figure(3);
    set(fighandle,'WindowButtonDownFcn',{@mytestcallback,x,y,z,M});
    %plot(f(50,1:N),'Color','red');
    %hold on;
    %plot(z(50,1:N),'Color','blue');
    imagesc(x,y,z);
    caxis([-1.5 1.5]);
    colorbar;
    
    %figure(1);
    %plot(z(50,1:N));
    %figure(2);
    %plot(z(50,1:N));    
end

function mytestcallback(hObject,eventdata,x,y,f,M)
    CP = get(get(hObject,'CurrentAxes'), 'CurrentPoint');
    CP(1,2);
    [~,xposition] = min(abs(x - CP(1,1)));
    figure(5);
    plot(y,f(1:M,xposition));
    
    hold on;
    ele = 1.60217646E-19;
    kB = 1.3806488E-23;
    
    A=1e17*(ele.^2)./(8.*kB)
    C=0.02*ele/(2*kB);
    
    %p = [0.5, 4, CP(1,2)];
    
    %fitfunc = inline(sprintf('%e*(p(1)./p(2))*(cosh((%e./p(2))*(x-p(3)))).^(-2)', A,C),'p','x');
    
    %y = f(yposition,1:N);
    %[p,r,j] = nlinfit(y,f(1:M,xposition),fitfunc,p);
    
    %figure(6);
    %plot (y, fitfunc(p,y));
    %hold off;
    %p
    
    ff=f(1:M,xposition);
    [pks,locs] = findpeaks(ff,'minpeakheight',max(ff)/4,'minpeakdistance',5);
    
    str = '';
    for i = 1:size(locs,2)
        p1 = (i-1)*3+1;
        p2 = (i-1)*3+2;
        p3 = (i-1)*3+3;
        p(p1) = pks(i);
        p(p2) = 4;
        p(p3) = y(locs(i));
        tmp = sprintf('%e*(p(%i)./p(%i))*(cosh((%e./p(%i))*(x-p(%i)))).^(-2)+', A,p1,p2,C,p2,p3)
        str = strcat(str,tmp);
    end
    str = strcat(str,'0');
    fitfunc = inline(str,'p','x');
    [p,r,j] = nlinfit(y,f(1:M,xposition),fitfunc,p);
    plot (y, fitfunc(p,y));
    ylim('manual');
    ylim([0 1]);
    hold on;
    %p = [0.5, 4, CP(1,1)];
    %plot (y, fitfunc(p,y));
    
    %s = fitoptions('Method','NonlinearLeastSquares');
    %f = fittype('a*(x-b)^n','problem','n','options',s); 
    %A*(gamma/T)*cosh^-1[(C/T)*(x-x0)]
    %[c2,gof2] = fit(cdate,pop,f,'problem',2);
    %plot(c2,'m');
end

