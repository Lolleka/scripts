
function [map] = CrNMR(T, newA, sig, plot_map)
    plot_pop = 0;
    plot_moments = 0;
    plot_3d = 0;
    freq_step = 0.5;
    freq_max = 50;
     
    global A;
    global h_lines;
    global f_lines;
    A = newA;
    Gamma = 2.406;
    
    %carica momenti
    moments_GS = load('.\Cr8Cd_moments_GS.dat');
    moments_ES1 = load('.\Cr8Cd_moments_ES1.dat');
    moments_ES2 = load('.\Cr8Cd_moments_ES2.dat');
    %carica i livelli
    levels = load('.\Cr8Cd_levels.dat');
    %carica dati sperimentali
    exp_pts = load('.\Cr8Cd_exp_points.dat');
    
    %livelli
    levels_h = levels(:,1);
    levels_GS = levels(:,2);
    levels_ES1 = levels(:,3);
    levels_ES2 = levels(:,4);
    %Spin totale
    %ST_GS = levels(:,4);
    %ST_ES1 = levels(:,5);
    %ST_ES2 = levels(:,6);
    
    %momenti
    Site_GS = moments_GS(:,2:5);
    Site_ES1 = moments_ES1(:,2:5);
    Site_ES2 = moments_ES2(:,2:5);

    %Interpola i livelli
    h = moments_GS(:,1)';
    ilevels_GS = interp1 (levels_h,levels_GS,h);
    ilevels_ES1 = interp1 (levels_h,levels_ES1,h);
    ilevels_ES2 = interp1 (levels_h,levels_ES2,h);
    %iST_GS = interp1 (levels_h,ST_GS,h);
    %iST_ES1 = interp1 (levels_h,ST_ES1,h);
    %iST_ES2 = interp1 (levels_h,ST_ES2,h);

    %calcola le frequenze di larmor
    Line_GS = abs((Site_GS*A+cat(2,h',h',h',h'))*Gamma);
    Line_ES1 = abs((Site_ES1*A+cat(2,h',h',h',h'))*Gamma);
    Line_ES2 = abs((Site_ES2*A+cat(2,h',h',h',h'))*Gamma);
    %Line_GS = abs((Site_GS*A)*Gamma);
    %Line_ES1 = abs((Site_ES1*A)*Gamma);
    
    %Funzione di partizione
    Z = 1+exp(-(ilevels_ES1-ilevels_GS)/T)+exp(-(ilevels_ES2-ilevels_GS)/T);;
    
    %Popolazioni
    Pop_GS = 1./Z;
    Pop_ES1 = exp(-(ilevels_ES1-ilevels_GS)/T)./Z; 
    Pop_ES2 = exp(-(ilevels_ES2-ilevels_GS)/T)./Z; 
    
    %plotta le popolazioni
    if plot_pop
        figure(1);
        p_GS = plot(h,Pop_GS);
        hold on;
        p_ES1 = plot(h,Pop_ES1);
        p_ES2 = plot(h,Pop_ES2);
        hold off;
        set(p_ES1,'Color','red');
        set(p_ES1,'Color','blue');
        ylim([0 1]);  
    end
    %plotta i momenti
    if plot_moments
        figure(2);
        plot (h,Line_GS);
        %figure(3);
        hold on;
        plot (h,Line_ES1,'--');
        plot (h,Line_ES2, '-.-');
        xlabel('Field (Tesla)');
        ylabel('Frequency (MHz)');
        hold off;
    end
    
    f = 0:0.1:freq_max;
    map=[];

    for i = 1:size(h,2)
        spectrum_GS =zeros(size(f));
        spectrum_ES1 =zeros(size(f));
        spectrum_ES2 =zeros(size(f));
        for il = 1:4
            spectrum_GS = spectrum_GS + gaussmf(f, [sig Line_GS(i,il)]);
            spectrum_ES1 = spectrum_ES1 + gaussmf(f, [sig Line_ES1(i,il)]);
            spectrum_ES2 = spectrum_ES2 + gaussmf(f, [sig Line_ES2(i,il)]);
        end
        map = [map; spectrum_GS*Pop_GS(i) + spectrum_ES1*Pop_ES1(i) + spectrum_ES2*Pop_ES2(i)];
    end
    
    if plot_map
        figure(3);
        imagesc(h,f,map');
        %Iimg = get(I,'CData');
        colormap(gray);
        text(0.05,0.8,['A = ' num2str(A)],'FontSize',18,'Units', 'normalized','Color',[1 1 1]);
        hold on;
        xlabel('Field (Tesla)');
        ylabel('Frequency (MHz)');
        set(gca,'YDir','normal');
        %set(fh,'WindowButtonMotionFcn',{@doslice,h,y,m});
        p_exp_pts = plot (exp_pts(:,2),exp_pts(:,1),'o');
        set(p_exp_pts,'Color','red');    
        hold off;
    end
    h_lines = h;
    f_lines = f;
    
    if plot_3d
        figure(4);
        [X,Y] = meshgrid(h, f);
        size(X)
        size(Y)
        size(m)
        surf(X,Y,m');
    end    
end


function doslice(hObject,eventdata,x,y,f)
    global A;
    CP = get(get(hObject,'CurrentAxes'), 'CurrentPoint');
    %CP(1,2)
    %CP(1,1)
    %[~,xposition] = min(abs(x - CP(1,1)));
    [~,yposition] = min(abs(y - CP(1,2)));
   % CP(1,2)
    figure(5);
    plot(x,f(:,yposition));
    xlabel('Field (Tesla)');
    ylabel('Signal Intensity (a.u.)');    
    text(0.05,0.9,[num2str(y(yposition)) ' MHz'],'FontSize',18,'Units', 'normalized');
    text(0.05,0.8,['A = ' num2str(A)],'FontSize',18,'Units', 'normalized');
    %hold on;
end


