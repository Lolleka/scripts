function testderivata (root)    

    files = dir (strcat(root,'*'));
    filename = files(1).name;

    close all;
    
    matrix = load(filename);    % battezzo "matrix" il filename che carico
    [row,col]=size(matrix);     % recupera le dimensioni della matrice del file dati ribattezzato "matrix"
                                % e inserisce questi due numeri (righe e colonne) nelle variabili row e col
    x=matrix(1,2:col);
    y=matrix(2:row,1);          % recuperano rispettivamente la prima riga e la prima colonna del mio file dati
                                % tralasciando il primo elemento, che ho settato a 0

    z=matrix(2:row,2:col);      % recupera il resto del mio file dati e lo inserisce in una matrice di nome z
    figure(1);                  % seleziona una finestra per fare la prima figura
    imagesc(x,y,z);             % usa i vettori x ed y per stampare in figura la mappa contenuta in z, producendo un plot a colori come in  origin
    caxis([-0.0125 0.0125]);          % definisce il range della scala di colore   
    colorbar
    %xlim([-1 1])                % setta il range dell'asse x
    
    d=[];       % non è inutile, battezza una nuova matrice di nome d, vuota
    row=row-1;
    col=col-1;      % diminuisce di 1 il numero di colonne e righe per le operazioni successive, 
                    % perchè la prima colonna e la prima riga del file dati sono state utilizzate per i vettori x ed y
                    % e non rientrano tra i dati di cui calcolare la derivata

     % il ciclo for successivo elabora tutte le righe della matrice z e fa le seguenti 2 cose:
     % calcola la derivata parziale di z rispetto a x con diff(z(i,1:col))./diff(x)
     %e assegna il risultato a un vettore di nome u, il cui primo elemento è 0
     % poi concatena questo vettore u alla matrice d, in modo da crearla riga per riga
     % così ho costruito la matrice contenente le derivate parziali dei miei dati, sulle righe               

    for i=1:row
        u=[0 diff(z(i,1:col))./diff(x)];
        d=[d; u];
    end
    figure(2);          % seleziona un secondo slot per plottare la matrice d
    imagesc(x,y,d);     % plotta la matrice d, ovvero la conduttanza differenziale
    %colorbar
    caxis([-0.0 0.3]);      
    %xlim([-1 1])
    
    saveas(1,[radice '_plot.jpg']);     % salva la figura 1 in formato jpg, col nome del file cioè "radice" seguito da "_plot"
    saveas(2,[radice '_der.jpg']);      % salva la figura 2 in formato jpg, col nome del file cioè "radice" seguito da "_der"
    
end