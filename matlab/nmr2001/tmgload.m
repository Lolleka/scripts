%TMGLOAD TecMaG binary data LOADer
%Carica i dati da files binari in formato TecMag (.tnt) in una matrice 
%di Matlab. Puo` caricare un singolo file 1-D o 2-D, oppure piu` files 1-D
%contenuti in una sottodirectory. 
%
%
%MODO FILE SINGOLO: 
%
%  [FID,SEQUENCE,BASE_FILENAME] = tmgload(FILENAME,option1,...,option3);
%
%dove la stringa FILENAME e` il nome di un singolo file .tnt 1-D o 2-D.
%
%I segnali di fid o di eco sono restituiti nella matrice complessa FID, di 
%dimensione [NPoints1D x NPoints2D] (ciascuna colonna di FID corrisponde 
%ad un esperimento unidimensionale).
%
%La matrice SEQUENCE contiene i dati rilevanti dell'acquisizione e della 
%sequenza di impulsi. La dimensione di SEQUENCE e` [(4 + 2*Npulses) x 
%NPoints2D], dove Npulses e` il numero di impulsi nella sequenza; ciascuna 
%colonna corrisponde a una colonna di FID.
%Le colonne di SEQUENCE hanno la struttura:
%
%SEQUENCE(:,k) = [NPoints1D; DwellTime; ReferenceFrequency; NScans; 
%	         PulseLength#1; DelayLength#1; ... ; 
%		 PulseLength#n; DelayLength#n]
% (tutti i tempi sono in secondi).
%
%Rispetto alla sequenze di impulsi completa nel file TecMag, tutti gli eventi 
%precedenti il primo impulso o successivi alla acquisizione vengono ignorati,
%ed eventi consecutivi non contenenti impulsi vengono fusi in un unico tempo
%di ritardo. Ad esempio, se tra il secondo impulso e l'acquisizione si trovano
%gli eventi RxBlanking = 10u, Acq_Delay = 5u, ad essi corrisponde 
%SEQUENCE(8,k) = 15e-6.
%
%BASE_FILENAME e` il nome del file letto, privato della parte di nome di 
%directory. Questo argomento di uscita e` implementato per ragioni di 
%omogeneita` con il modo multi-file.
%
%
%MODO MULTI FILE: 
%
%  [FID,sequence,base_filenames] = tmgload(FILENAME_PATTERN,option1,...,option3)
%
%dove la stringa FILENAME_PATTERN e` del tipo '<directory_name>/<pattern>',
%e <pattern> contiene una o piu` wildcards '*' o '@'.
%Ogni file 1-D in <directory_name> il cui nome combacia con <pattern> viene 
%caricato in FID. 
%Il carattere '*' ha il significato ordinario, e sostituisce qualsiasi 
%sottostringa. '@' si comporta come '*' con la restrizione che la 
%sottostringa sostituita sia interpretabile come un numero. 
%Ad esempio: se la directory /home/pippo/dati contiene i files
%
%  t40MHz.tnt
%  t_40MHz.tnt
%  u41MHz.tnt
%
%il comando
%
%  tmgload('/home/pippo/dati/t*.tnt') 
%
%carica t40MHz.tnt e t_40MHz.tnt, mentre
%
%  tmgload('/home/pippo/dati/t@.tnt')
%
%carica solo t40MHz.tnt perche' per 't_40MHz.tnt' la stringa sostituita '_40'
%non e` interpretabile come un numero. L'ordine di caricamento dei files e` 
%quello del dir listing, cioe` casuale.
%Di default, il confronto dei nomi e` case sensitive. Sotto Windows si puo` 
%rendere il confronto case insensitive impostando la variabile globale 
%TMGcasinsens ad un valore non nullo:
%  global TMGcasinsens;
%  TMGcasinsens = 1;
%La dimensione della matrice FID e` [Npt_max x Nfiles], dove Npt_max e` la 
%massima dimensione 1D dei files caricati, e Nfiles e` il numero dei files;
%ogni colonna di FID corrisponde ad un file. Se la dimensione 1D di un file
%e` inferiore a Npt_max, la colonna corrispondente viene terminata con zeri.
%
%Altri Argomenti di Output
%
%SEQUENCE ha lo stesso significato che nel modo singolo file. La dimensione
%e` [(4 + 2*Npulses_max) x Nfiles]. Se in qualche esperimento Npulses < 
%Npulses_max, la colonna corrispondente e` terminata con numeri negativi.
%
%BASE_FILENAMES (stringa) contiene i nomi dei files separati da 'a capo',
%nell'ordine in cui sono stati letti (l'm-esimo nome corrisponde alla
%m-esima colonna di FID e SEQUENCE).
%
%		________________________________________
%
%ARGOMENTI OPZIONALI (in entrambi i modi).
%Gli input di tipo stringa che iniziano con un '-' sono interpretati come 
%opzioni. Le opzioni possono essere inviate in qualsiasi ordine.
%Le opzioni valide sono '-n<scans>', '-t', '-v' .  
%
%La stringa '-n<scans>' attiva la normalizzazione dei FID a <scans> scansioni, 
%dove <scans> e` un intero. Ad esempio, se gli spettri sono stati ottenuti 
%mediando su 5000 scansioni e i dati vengono letti con l'opzione '-n100', i 
%fids vengono moltiplicati per 100/5000.
%L'opzione '-n' equivale a '-n1'.
%
%
%L'opzione '-t' attiva l'antitrasformata di Fourier automatica per i files 
%salvati dopo aver eseguito la FFT sul buffer dei dati. L'esecuzione di una o
%piu` trasformazioni e` segnalata da un flag nel file .TNT. Senza l'opzione 
%'-t' appare un messaggio di avviso ogni qual volta viene trovato alto il 
%flag, ma i dati non vengono antitrasformati.
%NB: 
%Il e` flag basso se e solo se NESSUNA trasformazione e` stata eseguita sul 
%file: in questo caso l'opzione '-t' e` sicura (l'antitrasformata non viene 
%eseguita).
%Il flag alto non garantisce pero` che i fid siano nel dominio delle frequenze.
%Ad esempio: 
% 1) nei files bidimensionali in cui solo alcune zone sono state trasformate;
% 2) nei files trasformati e antitrasformati prima di essere salvati;
%il flag rimane impostato. 
%Per distinguere dominio del tempo da dominio delle frequenze viene percio` 
%applicato il seguente criterio euristico (non affidabile al 100%): 
%i fid vengono ciclicamente antitrasformati fino a che l'ultimo canale del fid 
%e` un numero piccolo (compatibile con 0 in floating point a singola 
%precisione). Se il test fallisce dopo 4 trasformazioni (la FFT e` ciclica di 
%periodo 4), il programma da` un messaggio di avviso e il fid originale viene 
%comunque caricato.
% 
%
%L'opzione '-v' fa ritornare in FID il numero di versione del programma
%nella forma [major,minor, patchlevel]. Il programma non esegue nulla,
%e tutti gli altri argomenti vengono ignorati.
%
%		________________________________________
%
%LIMITAZIONI.
%1) Allo stato attuale TMGLOAD non supporta acquisizioni e tabelle di 
%   dimensione > 2
%2) Nella lettura delle sequenze, sono ignorati:
% -a) i loop;
% -b) gli offset di frequenza;
%
%					G.A. 27-10-00
%					ultimo aggiornamento 06-12-00
%					versione corrente: 1.4.2
%
