% template per funzione con numero variabile di argomenti (di input)

function comando(varargin)
    % Parte generale... come recuperare il numero dei parametri e i
    % parametri stessi... (FUNZIA??)
    files = dir (strcat(varargin{1},'*'));
    filename = files(1).name;
    fid = fopen(filename,'r');
    for nl = 1:10
        textline = fgetl(fid);
        fprintf('%s\n',textline);
    end
    fclose(fid);
end