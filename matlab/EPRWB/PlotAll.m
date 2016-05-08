legendlist = ['1'; '3'; '5'; '7'; '9'];
colorlist = prism(size(legendlist,1));
filelist_res = [legendlist  repmat(['spire_res.dat'],5,1)];
filelist_map = [legendlist  repmat(['spire_map.dat'],5,1)];
%figure(2);

%for i=1:size(filelist,1);
%    dat = load(filelist(i,:));
%    plot (dat(:,1),dat(:,2), 'Color', colorlist(i,:));
%    hold on;
%end

%legend(legendlist );
%hold off;

for i=1:size(filelist_res,1);
    fh = figure(2);
    clf;
    dat = load(filelist_res(i,:));
    [max_dat cI] = max(dat(:,2));
    plot (dat(:,1),dat(:,2)./max_dat, 'Color', 'red','LineWidth',2);
    hold on;
    plot (dat(:,1),dat(:,3), 'Color', 'blue');
    hold off;
    set(gca,'FontName','Arial','FontSize',14);
    %set(gca,'plotboxaspectratio',[4 1 1]);
    
    Coord=get(fh,'Position');
    NewLength=Coord(3)*4;
    set(fh,'Position',[Coord(1),Coord(2),1200,700]);

    set(get(gca,'XLabel'),'String','\fontname{Arial}\fontsize{14}Frequency(GHz)');
    set(get(gca,'YLabel'),'String','\fontname{Arial}\fontsize{14}Amplitude & Abs %');
    %title([legendlist(i) ' spire'],'FontName','Arial','FontSize',14);
    grid on;
    set(gcf,'PaperPositionMode','auto');
    legend('Amplitude','Abs %');
    %print('-dpsc2','-zbuffer','-r200');
    %print(fh,'-dpdf',[legendlist(i) ' spire.pdf'],'-zbuffer','-r200');
    print(fh,'-dpng',[legendlist(i) ' spire.png'],'-zbuffer','-r200');
    %legend('Amplitude', 'Abs/Disp Ratio');

    if 1;
    
    fh = figure(2);
    clf;
    datm = load(filelist_map(i,:));
    MCh1 = datm(2:end,2:end)';
    imagesc(datm(2:end,1),datm(1,2:end),MCh1);
    
    max_surf = max(max(MCh1));
    min_surf = min(min(MCh1));
    tot_span = max_surf-min_surf;
    n1 = floor(200*max_surf/tot_span);
    n2 = floor(200*abs(min_surf)/tot_span);
    mymap = [];
    mymap =  [mymap; [linspace(0,1,n2)' linspace(0,1,n2)' linspace(1,1,n2)']];
    mymap =  [mymap; [linspace(1,1,30)' linspace(1,1,30)' linspace(1,1,30)']];
    mymap =  [mymap; [linspace(1,1,n1)' linspace(1,0,n1)' linspace(1,0,n1)']];
    colormap(mymap);
    %caxis([-1.5 1.5]);
    colorbar;  
    
    ax1 = gca;
    set(gca,'FontName','Arial','FontSize',14);
    set(get(gca,'XLabel'),'String','\fontname{Arial}\fontsize{14}Frequency(GHz)');
    set(get(gca,'YLabel'),'String','\fontname{Arial}\fontsize{14}Field (Tesla)');
    %title([legendlist(i) ' spire'],'FontName','Arial','FontSize',14);    
    %hold off;
    %figure(3);
    ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
       
    set(ax2,'xtick',[],'ytick',[]);
    hold on;
    
    plot(dat(:,1),dat(:,3), 'color', 'blue','Parent',ax2);
    
    
    linkaxes([ax1 ax2],'x');
    xlim(ax1,[2 20]);
    xlim(ax2,[2 20]);
    ylim(ax2,[0 5]);
    set(fh,'CurrentAxes',ax2);
    line([2 20],[1 1],'LineStyle','--','Color','black');
    
    set(gcf,'PaperPositionMode','auto');
    %print('-dpsc2','-zbuffer','-r200');
    %print(fh,'-dpdf',[legendlist(i) ' spire_map.pdf'],'-zbuffer','-r200');    
    print(fh,'-dpng',[legendlist(i) ' spire_map.png'],'-zbuffer','-r200');    
    
    hold off;
    %colorbar('delete');
    end
    %title ('Ch1 intensity map - Frequency vs Field');
    %hold on;
    %plot (dat(:,1),dat(:,3), 'Color', 'blue');
    %hold off;
    %set(gca,'FontName','Arial','FontSize',14);
    %set(gca,'plotboxaspectratio',[4 1 1]);
    
    %Coord=get(fh,'Position');
    %NewLength=Coord(3)*4;
    %set(fh,'Position',[Coord(1),Coord(2),1200,700]);
    
    %set(get(gca,'XLabel'),'String','\fontname{Arial}\fontsize{14}Frequency(GHz)');
    %set(get(gca,'YLabel'),'String','\fontname{Arial}\fontsize{14}Amplitude & Abs/Disp Ratio');
    %title([legendlist(i) ' spire'],'FontName','Arial','FontSize',14);
    %grid on;
    %set(gcf,'PaperPositionMode','auto');
    %print('-dpsc2','-zbuffer','-r200');
    %print(fh,'-dpng',[legendlist(i) ' spire.png'],'-zbuffer','-r200');
    %legend('Amplitude', 'Abs/Disp Ratio');
end
