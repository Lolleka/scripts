function BAScript(base_path)

    f = dir([base_path '*.csv']);
    Ladder_i = find(1-cellfun('isempty',regexp({f.name},'\w*Ladder*','match','once')) == 1);
    Ladder_name = f(Ladder_i).name;
    f(Ladder_i) = [];
    
    %Allocate an empty matrix variable to stare Curve_list data
    Curve_list={};
    cmap = hsv(length(f));
    for i = 1:length(f)
        Curve_list{i} = {BALoad([base_path f(i).name],[base_path Ladder_name]) f(i).name cmap(i,:)};
    end

    %Create new figure
    hfig = figure(1);
    %Aux variable to store curve y maximum
    max_y = 0;
    %Loop over curve list
    for i = 1:length(Curve_list)
        %Get i-th curve
        Curve = Curve_list{i}{1};
        %Renormalize Curve, i.e., divide by total area
        Curve(:,2) = Curve(:,2)./trapz(Curve(:,1),Curve(:,2));
        %Update maximum (max_y)
        tmp_max = max(Curve(:,2));
        if tmp_max > max_y
            max_y = tmp_max;
        end
        %Plot curve, get color from cell matrix (i-th row, 3rd column)
        plot(Curve(:,1),Curve(:,2),'Color',Curve_list{i}{3});
        save([base_path Curve_list{i}{2} '.txt'], 'Curve', '-ascii');
        %Do not clear graph on next plot
        hold on;    
    end
    %Allow cleanup on next plot
    hold off;

    XScale_value = 'lin';
    %sets the current XScale to either linear (lin) or logarithmic (log)
    %sets the Y range (0 to max_y*1.2), for some room
    %sets the X range to fit the curve to the plot window
    set(gca,'XScale',XScale_value,'YLim',[0 max_y*1.2],'XLim',[min(Curve_list{1}{1}(:,1)) max(Curve_list{1}{1}(:,1))]);
    %sets x and y labels
    xlabel('Base Pairs');
    ylabel('\rho (BP)');

    %Get Legend tags
    tmp = reshape([Curve_list{:}],3,length(Curve_list));
    %Create Legend
    hleg = legend(tmp{2,:});
end