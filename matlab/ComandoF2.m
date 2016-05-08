% template per funzione con numero variabile di argomenti (di input)

function comando(varargin)

    % Parte generale... come recuperare il numero dei parametri e i
    % parametri stessi... (FUNZIA??)
    fprintf('Chiamata a funzione con %d parametri.\n',nargin);
    fprintf('Lista parametri:\n');
    for i=1:nargin
        fprintf('#%d - "%s"\n',i,varargin{i});
    end
    fprintf('\n');
    
    % nel caso specifico, assummo che il parametro (se esiste) sia il file
    % delle misure. Mi aspetto qualcosa del tipo m###.dat
    
    pathname = [pwd '\'];     % retreive current pathname
    if nargin>0
        filename = varargin{1};
        if filename(max(1,end-3))~='.'
            if filename(1)=='m'
                % case m### => append .dat
                filename = [filename '.dat'];
            else
                % just a number => recreate filename
                filenumber = str2num(filename);
                filename   = sprintf('m%3d.dat',filenumber);
                spaces = find(filename==' ');
                filename(spaces)='0';
            end
        end 
        %  check if the file really exists!
        if exist(filename,'file')==0
            [filename pathname] = uigetfile('*.dat','File non trovato: scegli il file da caricare');
        end
    else
        % nessun parametro, apri GUI per scegliere il file
        [filename pathname] = uigetfile('*.dat','Scegli il file da caricare');
    end
    
    fprintf('Path: %s\nFile: %s\n\nSeguono prime righe del file...\n\n',pathname,filename);
    
    % Qui continua l'operazione... ora apre il file e stampa le prime dieci
    % righe.
    
    fid = fopen(filename,'r');
    for nl = 1:10
        textline = fgetl(fid);
        fprintf('%s\n',textline);
    end
    fclose(fid);
    
end