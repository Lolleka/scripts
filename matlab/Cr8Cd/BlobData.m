
%Variabili globali per la mappa teorica in arrivo da CrNMR
global lines_map;
global h_lines;
global f_lines;

lines_map = CrNMR(1.35,-10.8, 1,0);

global bdata;

Nh = 200;
Nf = 150;
hmax = 8;
fmax = 40;
h = 0:(hmax/Nh):hmax;
f = 0:(fmax/Nf):fmax;

sweeps = {'Hsweep_16MHz.dat' 'Hsweep_20MHz.dat' 'Hsweep_23MHz.dat' 'Hsweep_25MHz.dat' 'Hsweep_32MHz.dat' 'Fsweep_1.8T.dat' 'Fsweep_5T.dat' 'Fsweep_5.7T.dat' 'Fsweep_8T.dat'};
%data_map = mat2gray(MapData(sweeps{1}, h , f));
data_map = MapData(sweeps{1}, h , f);

for i = 2:length(sweeps)
    %data_map = data_map + mat2gray(MapData(sweeps{i}, h , f));
    data_map = data_map + MapData(sweeps{i}, h , f);
end
data_map = mat2gray(data_map);

[X,Y] = meshgrid(f, h);  
result = [];
%TotNsig = 15;
TotNA = 9.5;
Astart = 40;
Aend = 11.5;
Astep = (Aend - Astart)/(TotNA-1);
%sigstart = 0.01;
%sigend = 0.1;
%sigstep = (sigend - sigstart)/(TotNsig-1);
sig = 2.1
A = -10.6

do_optimization = 0

if do_optimization
    wb = waitbar(0, 'Finding optimal coupling constant and line-width...')

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
%figure(2);
%RX = Astart:Astep:Aend;
%RY = sigstart:sigstep:sigend;
%RM = reshape(result(:,1),length(RY),length(RX));
%surf(RX,RY,RM);

%A
%sig
lines_map = CrNMR(1.35,A,sig,0);
%lines_map = CrNMR(1.35,-10.5,1.8,0);

fig_h = figure(1);
imagesc(h_lines,f_lines,lines_map');
colormap(gray);

green = cat(3, zeros(size(data_map')),ones(size(data_map')), zeros(size(data_map')));
%hold on;
%green_overlay = imagesc(h,f,green); 
%hold off;
%set(green_overlay, 'AlphaData', mat2gray(data_map'));
xlabel('Field (Tesla)');
ylabel('Frequency (MHz)');
set(gca,'YDir','normal');
grid;
text(0.05,0.8,['A = ' num2str(A) '; sig =' num2str(sig)],'FontSize',18,'Units', 'normalized','Color',[1 1 1]);

saveas (fig_h, '.\pic_opti.png', 'png');