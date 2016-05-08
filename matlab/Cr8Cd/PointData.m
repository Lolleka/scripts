function PointData(A)
%Variabili globali per la mappa teorica in arrivo da CrNMR
global lines_map;
global h_lines;
global f_lines;

%lines_map = CrNMR(1.35,-10.8, 1,0);

global bdata;

fid = fopen('exp_points.dat','rt');
%Skip header and parse file
datacell = textscan(fid, '%f%f%f%f%f', 'Delimiter', ',', 'HeaderLines', 1, 'CollectOutput', 1);
%Close file
fclose(fid);
%Return data
exp_points = datacell{1};
exp_points = exp_points(find(exp_points(:,5)==1),:);
%Cleanup
clear datacell;

exp_limits = load('exp_limits.dat');
%exp_points = PeakData('Fsweep_1.8T.dat', 'f')
%exp_points  = [exp_points; PeakData('Fsweep_5T.dat', 'f')];
%exp_points  = [exp_points; PeakData('Fsweep_5.7T.dat', 'f')];
%exp_points  = [exp_points; PeakData('Hsweep_20Mhz.dat', 'h')];
%exp_points  = [exp_points; PeakData('Hsweep_25Mhz.dat', 'h')];
%result = [];
%TotNsig = 15;
%TotNA = 9.5;
%Astart = 40;
%Aend = 11.5;
%Astep = (Aend - Astart)/(TotNA-1);
%sigstart = 0.01;
%sigend = 0.1;
%sigstep = (sigend - sigstart)/(TotNsig-1);


do_optimization = 0;

figure(1);
mat_lines=CrLines(A);

%exp_points = exp_points(1:7,:);
%exp_points = exp_points(8:17,:);

for i = 1:size(exp_limits,1)
    line(exp_limits(i,1:2), exp_limits(i,3:4),[0 0], 'Color', [0.5 0.5 0.5]);
end

hold on;
plot(exp_points(:,1),exp_points(:,2),'o');
hold on;
for i = 1:size(exp_points,1)
    pt = exp_points(i,:);
    X = [pt(1)-pt(3) pt(1)+pt(3)];
    Y = [pt(2) pt(2)];
    line (X,Y,[0 0]); %,'MarkerSize',4);
    X = [pt(1) pt(1)];
    Y = [pt(2)-pt(4) pt(2)+pt(4)];
    line (X,Y);
end

p = plot(mat_lines(:,1),mat_lines(:,2:5),'-');
set(p,'color','black');
p = plot(mat_lines(:,6),mat_lines(:,7:10),'-');
set(p,'color','red');
p = plot(mat_lines(:,11),mat_lines(:,12:15),'-');
set(p,'color','blue');
p = plot(mat_lines(:,16),mat_lines(:,17:20),'-');
set(p,'color','green');
hold off;
set(gca, 'YLim', [18 38], 'XLim', [0 9]);
text(0.05,0.8,['A = ' num2str(A)],'FontSize',18,'Units', 'normalized','Color','black');
%mat_lines_GS = mat_lines(:,1:6);
%mat_lines_ES1 = [mat_lines(:,1) mat_lines(:,7:11)];
%mat_lines_ES2 = [mat_lines(:,1) mat_lines(:,12:16)];

%mat_lines = [mat_lines_GS; mat_lines_ES1; mat_lines_ES2];

%save('exp_points.dat','exp_points','-ascii');
save('mat_lines.dat','mat_lines','-ascii');
%c=0;
%pc = 0;
%mat = [[][][]];
%for i = 1:size(mat_lines,1)
%    if mat_lines(i,6) == 0
%        mat(1) = [mat(1); mat_lines(i,:)];
%        c=0;
%    elseif mat_lines(i,6) == 1
%        mat(2) = [mat(2); mat_lines(i,:)];
%        c=1;
%    elseif mat_lines(i,6) == 2
%        mat(3) = [mat(3); mat_lines(i,:)];
%        c=2;
%    end
%    
%    if pc ~= c
%        mat(pc) = [mat(pc); repmat(0,1,6)];
 %       pc = c;
 %   end
%end

%hold on;

%for i = 2:size(mat_lines,2)
%    ind = find(diff(diff(mat_lines(:,i)')) > 0.1);
%    mat_lines(ind+2,:) = mat_lines(ind+1,:);
%    mat_lines(find(diff(diff(mat_lines(:,i)')) > 0.1)+1,i) = nan;
%    p = plot(mat_lines(:,1),mat_lines(:,i),'-');
%end

%hold off;

if do_optimization
    wb = waitbar(0, 'Finding optimal coupling constant and line-width...');
    for i= 1:TotNA
        %for j = 1:TotNsig
            A = Astart+Astep*i;
            %sig =sigstart+sigstep*j;
            %P = (j + (i-1)*TotNsig)/(TotNsig*TotNA);
            P = i/TotNA;
            wb = waitbar(P)
            lines_map = CrNMR(1.35, -A, sig, 0);
            z = interp2(f_lines,h_lines,lines_map,X,Y);
            mu = (z - data_map).^2;
            result = [result; sum(mu(:)) A sig];
        %end
    end
    close (wb)
    opti = result(find(result(:,1) == max(result(:,1))),:);
    A = -opti(2);
    sig = opti(3);
end

%lines_map = CrNMR(1.35,A,sig,0);
%lines_map = CrNMR(1.35,-10.5,1.8,0);
%saveas (fig_h, '.\pic_opti.png', 'png');
end