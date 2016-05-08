
[num,txt,raw] = xlsread('lista.xlsx','materiale');
nans = cellfun(@(V) any(isnan(V(:))), raw);
raw(find(nans(:,1) == 1),:) = [];
nans = cellfun(@(V) any(isnan(V(:))), raw);
raw(find(nans(:,3) == 1),3) = {[0]};

header = raw(1,:);
raw(1,:) = [];

%now raw contains only useful data

Group=raw(:,1);
Item=raw(:,2);
Weight=cell2mat(raw(:,3));
P = {};
ps = find(nans(1,5:end) == 0)+4;


nans = cellfun(@(V) any(isnan(V(:))), raw);
raw(find(nans(:,4) == 1),4) = {[0]};
Common=cell2mat(raw(:,4));
for i = 1:length(ps)
    raw(find(nans(:,ps(i)) == 1),ps(i)) = {[0]};
    P{1,i} = header{ps(i)};
    P{2,i}=cell2mat(raw(:,ps(i)));
end

tot = sum(Weight'.*(Common'+sum(cell2mat(P(2,:))')));
persone = size(P,2);
partot=[];

for i = 1:persone
    partot(i) = sum(Weight.*P{2,i});
end
partot


[S,ind] = sortrows(Weight);
ind=flipud(ind);
Weight = Weight(ind);
Item = Item(ind);
Group = Group(ind);
Common = Common(ind);
for i = 1:persone
    P{2,i} = P{2,i}(ind);
end

%for i = 1:length(Weight)
%    if Selected(i)
%        fprintf('%s\t%s\t%i\n',Group(i,:),Item(i,:),Weight(i))
%    end
%end
%fprintf('...........................\n')
%fprintf('Totale\t\t\t\t%i\n',tot)
%persone = 2;
%p=cell(length(P),1);

for i = 1:persone
    P{3,i} = [];
end

s=partot;

for i = 1:length(Weight)
    for p = 1:persone
        for j = 1:P{2,p}(i)
            P{3,p} = [P{3,p} i];
        end
    end
    %fprintf('%s\n',cell2mat(Item(i)));
    for j = 1:Common(i)
        [m,l] = min(s);
        P{3,l} = [P{3,l} i];
        s(l) = s(l)+Weight(i);
        
    end    
end


%Crea le liste:

%Legenda cell array P
% P{1,:} Riga nomi persone
% P{2,:} contiene gli array che definiscono la lista degli item propri di
% una persona
% P{3,:} contiene gli array che definiscono le liste elaborate
% dall'algoritmo (con ripetizioni in un primo momento, poi compressi prima
% di creare la lista. Ovvero P{3,:}(1,:) contiene gli elementi unici e P{3,:}(2,:)
% contiene i moltiplicatori relativi alla riga precedente)
% P{4,:} contiene gli array delle classi relative agli oggetti in P{3,:}
%Lista per persona:

fID = fopen('Lista.tex','w');

fprintf(fID, '\\documentclass[]{article}\n');

fprintf(fID, '\\usepackage{amsmath}\n');
fprintf(fID, '\\usepackage{amssymb}\n');
fprintf(fID, '\\usepackage{array}\n');
fprintf(fID, '\\usepackage{graphicx}\n');
fprintf(fID, '\\usepackage{url}\n');
fprintf(fID, '\\usepackage{multirow}\n');
fprintf(fID, '\\usepackage{epstopdf}\n');
fprintf(fID, '\\usepackage{bm}\n');
fprintf(fID, '\\usepackage{color}\n');
fprintf(fID, '\\usepackage{colortbl}\n');
fprintf(fID, '\\usepackage{xspace}\n');
fprintf(fID, '\\usepackage[italian]{babel}\n');
fprintf(fID, '\\usepackage[T1]{fontenc}\n');
fprintf(fID, '\\usepackage[utf8]{inputenc} \n');

%fprintf(fID, '\\definecolor{name}{system}{definition}\n');

fprintf(fID, '\\definecolor{Gray}{gray}{0.9}\n');
fprintf(fID, '\\definecolor{White}{gray}{1}\n');
fprintf(fID, '\\definecolor{Yellow}{rgb}{1,1,0.6}\n');
fprintf(fID, '\\definecolor{LightCyan}{rgb}{0.88,1,1}\n');


fprintf(fID, '\\title{Lista trekking}\n');
%fprintf(fID, '\\author{}\n');
fprintf(fID, '\\begin{document}\n');
fprintf(fID, '\\maketitle\n');

%fprintf(fID, '\\color{red}\n');

fprintf(fID, '\\centering\n');
fprintf(fID, '\\large\n');

fprintf(fID,'\\begin{tabular}{|l|r|}\n');
fprintf(fID,'\\hline\n');
fprintf(fID,'\\multicolumn{2}{|c|}{%s}\\\\\n', 'Panoramica pesi');
fprintf(fID,'\\hline\n');

flag = false;
for i=1:persone
    if flag
        col = 'Yellow';
    else
        col = 'White';
    end
    flag = ~flag;
    
    fprintf(fID,'\\rowcolor{%s}\n',col);
    fprintf(fID, 'Totale %s & %u g\\\\\n', P{1,i}, s(i));
    fprintf(fID,'\\rowcolor{%s}\n',col);
    fprintf(fID, 'Parziale solo %s & %u g\\\\\n', P{1,i}, partot(i));
    fprintf(fID,'\\rowcolor{%s}\n',col);
    fprintf(fID, 'Parziale comune %s & %u g\\\\\n', P{1,i}, s(i)-partot(i));
end

fprintf(fID,'\\hline\n');
fprintf(fID,'\\end{tabular}\n\n');

for i=1:persone
    %Ricalcola l'array lista nel formato item * moltiplicatore
    newP  = unique(P{3,i});
    
    if ~isempty(newP)
        countP = histc(P{3,i},newP);
        P{3,i} = [newP; countP];
        P{4,i} = Group(newP);
        [S, ind] = sortrows(P{4,i});
        P{4,i} = P{4,i}(ind);
        P{3,i} = P{3,i}(:,ind);
        fprintf(fID,'\\section{%s (Totale: %u g)}\n\n', P{1,i}, s(i));
        
        %Crea la lista completa
        PrevGroup = cell2mat(P{4,i}(1));
        fprintf(fID,'\\begin{tabular}{|c|p{9cm}>{\\raggedleft\\arraybackslash}p{1.5cm}c|}\n');
        fprintf(fID,'\\hline\n');
        fprintf(fID,'\\multicolumn{4}{|c|}{%s}\\\\\n', PrevGroup);
        fprintf(fID,'\\hline\n');
        subtot=0;
        tot=0;
        flag = false;
        for j = 1:size(P{3,i},2)
             if ~strcmpi(PrevGroup,cell2mat(P{4,i}(j)))
                PrevGroup = cell2mat(P{4,i}(j));
                fprintf(fID,'\\rowcolor{LightCyan}\n');
                fprintf(fID,' & Subtotale Categoria & %u g & \\\\\n', subtot);
                fprintf(fID,'\\hline\n');
                fprintf(fID,'\\end{tabular}\n\n');
                fprintf(fID,'\\begin{tabular}{|c|p{9cm}>{\\raggedleft\\arraybackslash}p{1.5cm}c|}\n');
                fprintf(fID,'\\hline\n');
                fprintf(fID,'\\multicolumn{4}{|c|}{%s}\\\\\n', PrevGroup);
                fprintf(fID,'\\hline\n');
                subtot=0;
             end

            multiplier = P{3,i}(2,j);
            itemname = cell2mat(Item(P{3,i}(1,j)));
            totalitemweight = Weight(P{3,i}(1,j))*P{3,i}(2,j);
            if flag
                fprintf(fID,'\\rowcolor{Yellow}\n');
            end
            flag = ~flag;
            
            fprintf(fID,'%ux & %s & %u g & $\\square$\\\\\n', multiplier, itemname, totalitemweight);
            subtot = subtot + totalitemweight;
            tot = tot + totalitemweight;
        end
        fprintf(fID,'\\rowcolor{LightCyan}\n');
        fprintf(fID,' & Subtotale Categoria & %u g & \\\\\n', subtot);
        fprintf(fID,'\\hline\n');
        fprintf(fID,'\\end{tabular}\n\n');
        fprintf(fID,'\\newpage\n');
    end    
end
fprintf(fID, '\\end{document}\n');

fclose(fID);

!cat Lista.tex
!pdflatex Lista.tex