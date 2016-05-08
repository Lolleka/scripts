function PlotDPPH(filename)
    M = load(filename);
    M = M(:,2:6);
    pos = [1; find(diff(M(:,1))>0)];
    range = [pos(1:(length(pos)-1)) pos(2:length(pos))];
    range = [range; [pos(end) size(M,1)]];
    Hrange = [-0.001:0.00001:0.001];
    %Hrange = [0:0.0001:0.02];
    MCh1 = [];
    MCh2 = [];
    f = [];
    r = [];
    
    for i = 1:size(range,1)
        f = [f; M(range(i,1)+1,1)];    
        H = M(range(i,1):range(i,2),2);
        Ch1 = M(range(i,1):range(i,2),3);
        Ch2 = M(range(i,1):range(i,2),4);

        %T = [H smooth(Ch1,0.01) smooth(Ch2,0.01)];
        T = [H Ch1 Ch2];
        [res,pos] = sort(H);
        T = T(pos,:);
        T(find(diff(T(:,1))==0),:) = [];

        T(1:20,:) = [];        
        
        %Trova il livello di rumore medio entro i primi 20 canali
        %noise = sum(abs(T(1:30,2)))/30;
        cumul = cumtrapz(T(:,1),abs(T(:,2)));
        tmpy = T(:,1);
        tmpy(find(diff(cumul)==0),:) = [];
        cumul(find(diff(cumul)==0),:) = [];
        
        cumul = cumul/max(cumul);
        %figure(5); plot(tmpy,cumul);
        %find(diff(T(:,1))==0)
        central = interp1(cumul,tmpy,0.5);
        [minCh1,Imin] = min(T(:,2));
        [maxCh1,Imax] = max(T(:,3));    

        T(:,1) = T(:,1) - central;
        %figure(4);
        %plot(T(:,1), T(:,2),'color','black');
        
        [rtmp p] = GetRatio(T(:,1),T(:,2));
        r = [r; rtmp];
        central = p(3);
        T(:,1) = T(:,1) - central;
        
        %plot(T(:,1), T(:,2),'color','black');
        
        T = [T(:,1) smooth(T(:,2),0.03) smooth(T(:,3),0.03)];

        newCh1=interp1(T(:,1),T(:,2),Hrange);
        newCh2=interp1(T(:,1),T(:,3),Hrange);
        MCh1 = [MCh1; newCh1];
        MCh2 = [MCh2; newCh2];
    end
    %[f c]
    %fighandle = figure(1);
    %plot(f(50,1:N),'Color','red');
    %hold on;
    %plot(z(50,1:N),'Color','blue');
    MCh1(isnan(MCh1)) = 0;
    MCh2(isnan(MCh2)) = 0;

    %rows = find(f > 14);
    %f(rows) = [];
    %MCh1(rows,:) = [];

    %imagesc(Hrange,f,MCh1);

    %caxis([-5e-6 5e-6]);

    %colorbar;
    
    f_n = f(1):(f(end)-f(1))/400:f(end);
    Hrange_n = Hrange(1):(Hrange(end)-Hrange(1))/400:Hrange(end);
    [xi, yi] = meshgrid(Hrange_n,f_n);
    
    Hrange_n = Hrange;
    f_n = f;
    %f(1):0.01:f(end)
    %MCh1 = interp2(Hrange, f, MCh1, xi, yi, 'spline');

    %camlight right
    %set(gca, 'CameraPosition', [45 35 9.8])
    %mymap =  [linspace(1,1,100)' linspace(0,0,100)' linspace(0,0,100)'];
    max_surf = max(max(MCh1));
    min_surf = min(min(MCh1));
    tot_span = max_surf-min_surf;
    n1 = floor(200*max_surf/tot_span);
    n2 = floor(200*abs(min_surf)/tot_span);
    mymap = [];
    mymap =  [mymap; [linspace(0,1,n2)' linspace(0,1,n2)' linspace(1,1,n2)']];
    mymap =  [mymap; [linspace(1,1,30)' linspace(1,1,30)' linspace(1,1,30)']];
    mymap =  [mymap; [linspace(1,1,n1)' linspace(1,0,n1)' linspace(1,0,n1)']];

    %mymap = [];
    %mymap =  [mymap; [linspace(0,0,n2)' linspace(0,0,n2)' linspace(1,0,n2)']];
    %mymap =  [mymap; [linspace(0,1,n1)' linspace(0,0,n1)' linspace(0,0,n1)']];

    %mymap = [];
    %mymap =  [mymap; [linspace(0,0,n2)' linspace(0,0,n2)' linspace(1,0,n2)']];
    %mymap =  [mymap repmat([1 1 1], n2, 1)];

    %mymap = [ [1 0 0]; [0 0 1] ];
    %mymap = repmat (mymap,100,1);
    fighandle = figure(1); %,'Renderer','opengl');
    
    set(fighandle,'WindowButtonMotionFcn',{@MouseCallback,f_n,Hrange_n,MCh1',[min_surf max_surf]*1.1});
    imagesc(f_n,Hrange_n,MCh1');
    %caxis([-1.5 1.5]);
    colorbar;    
    title ('Ch1 intensity map - Frequency vs Field');
    %surf(Hrange_n,f_n,MCh1,'EdgeColor','none'); %,'MeshStyle','row');
    %shading interp;

    %colormap(mymap);
    %alpha(.8)
    %colormap([hsv(200);flipud(hsv(200))]);
    %colormap(hsv(200));
    
    a = [];
    for i = 1:length(f_n)
        %MCh1(i,:)
        a = [a; trapz(Hrange_n,abs(MCh1(i,:)))];
    end

    fi = figure(2);
    clf;
    %hold off;
    plot(f_n, a, 'color', 'red','LineWidth',3);
    title (['Ch1 intensity profile - Intensity vs Frequency - ' filename]);
    ax1 = gca;
    %hold off;
    %figure(3);
    ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
       
    set(ax2,'xtick',[]);
    hold on;
    plot(f_n, r, 'color', 'blue','Parent',ax2);
    
    
    linkaxes([ax1 ax2],'x');
    hold off;
    %xlim([2 9]);
    print(fi,'-djpeg',[filename '.jpg']);
    %title ('Ch1 intensity profile - Abs/Disp Ratio vs Frequency');
    
    %savedata = sortrows(savedata,1);
    savedata = [f_n a r];
    save(strrep(filename,'.dat','_res.dat'),'savedata','-ascii');
    savedata = [ [0 Hrange_n]; [f_n MCh1]];
    
    save(strrep(filename,'.dat','_map.dat'),'savedata','-ascii');
end

function MouseCallback(hObject,eventdata,x,y,z,ylims)
    CP = get(get(hObject,'CurrentAxes'), 'CurrentPoint');
    CP(1,2);
    [~,xposition] = min(abs(x - CP(1,1)));
    
    p = FitEPR(y',z(:,xposition));
    func_abs = inline ('-p(1).*(p(2).*(x-p(3)))./(pi.*(p(2).^2/4+(x-p(3)).^2).^2)','p','x');
    func_disp = inline ('-p(4).*(16.*p(2).*(p(2).^2-12.*x.^2-12.*p(3).^2+24.*x.*p(3)))./(pi.*(p(2).^2+4*x.^2+4.*p(3).^2-8.*x.*p(3)).^3)','p','x');
    funcstr =  '-p(1).*(p(2).*(x-p(3)))./(pi.*(p(2).^2/4+(x-p(3)).^2).^2)-p(4).*(16.*p(2).*(p(2).^2-12.*x.^2-12.*p(3).^2+24.*x.*p(3)))./(pi.*(p(2).^2+4*x.^2+4.*p(3).^2-8.*x.*p(3)).^3)';
    fitfunc = inline(funcstr,'p','x');    
    
    %GetRatio(y',z(:,xposition));
    
    figure(4);
    plot(y,z(:,xposition),'color','black');
    
    
    hold on;
    plot(y,fitfunc(p,y),'color','red');
    plot(y,func_abs(p,y),'color','green');
    plot(y,func_disp(p,y),'color','blue');
    t = title (['Signal at ' num2str(x(xposition),'% 2.2f') ' GHz'],'FontSize',20);
    set(t, 'FontSize', 26);
    %ylim (ylims);
    %set(gca, 'YLim'
    hold off;
   
end

function [ratio best_p] = GetRatio(x,f)
    p = FitEPR(x,f);
    func_abs = inline ('-p(1).*(p(2).*(x-p(3)))./(pi.*(p(2).^2/4+(x-p(3)).^2).^2)','p','x');
    func_disp = inline ('-p(4).*(16.*p(2).*(p(2).^2-12.*x.^2-12.*p(3).^2+24.*x.*p(3)))./(pi.*(p(2).^2+4*x.^2+4.*p(3).^2-8.*x.*p(3)).^3)','p','x');
    funcstr =  '-p(1).*(p(2).*(x-p(3)))./(pi.*(p(2).^2/4+(x-p(3)).^2).^2)-p(4).*(16.*p(2).*(p(2).^2-12.*x.^2-12.*p(3).^2+24.*x.*p(3)))./(pi.*(p(2).^2+4*x.^2+4.*p(3).^2-8.*x.*p(3)).^3)';
    fitfunc = inline(funcstr,'p','x');
    
    xx = x(1):((x(end)-x(1))/1000):x(end);
    
    %interpolate disp_diff 
    disp_diff_i = interp1(x,func_disp(p,x),xx);
    %integrate disp
    disp_cum = cumtrapz(xx, disp_diff_i);
    %hilbert transform to get disp
    abs_h = imag(hilbert(disp_cum));
    
    abs_h_diff = diff(abs_h)./diff(xx);
    abs_h_diff = [abs_h_diff 0]; 
    %abs_int = trapz(x,abs(func_abs(p,x)));
    %disp_int = trapz(x,abs(func_disp(p,x)));
    %ratio = abs_int/(abs_int+disp_int);    
    
    
    %figure(4);
    %hold off;
    
    %plot(y,z(:,xposition),'color','black');
    %plot(x,f,'color','black');
    %plot(x,func_abs(p,x),'color','red');
    
    %hold on;
    
    %plot(x,func_disp(p,x),'color','green');
    %plot(xx,abs_h_diff,'color','blue');
    %plot(x,func_abs(p,x),'color','green');
    
    
    abs_int = trapz(x,abs(func_abs(p,x)));
    abs_h_int = trapz(xx,abs(abs_h_diff));
    ratio = abs_int /(abs_h_int + abs_int) ;
    best_p = p;
end

function [best_p] = FitEPR(x,f)
    max_f = max(f);
    %min_f = abs(min(f));
    
    p = [2E-11 0.001 0 4E-15];
    p(1) = p(1) * max_f/ 2E-5;
    p(4) = p(4) * max_f/ 2E-5;
    
    p1 = p;
    p2 = p;
    p2(1) = -p2(1);
    p2(4) = -p2(4);
    
    funcstr =  '-p(1).*(p(2).*(x-p(3)))./(pi.*(p(2).^2/4+(x-p(3)).^2).^2)-p(4).*(16.*p(2).*(p(2).^2-12.*x.^2-12.*p(3).^2+24.*x.*p(3)))./(pi.*(p(2).^2+4*x.^2+4.*p(3).^2-8.*x.*p(3)).^3)';
    fitfunc = inline(funcstr,'p','x');
    [p1,r1,j] = nlinfit(x,f,fitfunc,p1);
    %p(4) = -p(4)
    %p(1) = p(1)*1.2;
    [p2,r2,j] = nlinfit(x,f,fitfunc,p2);
    
    chi1 = sum(r1.^2);
    chi2 = sum(r2.^2);
    
    if chi1 < chi2
        best_p = p1;
    else
        best_p = p2;
    end
end
