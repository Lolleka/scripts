function Digitize(path)
    fh = figure (1);
    %set (fh,'units','normalized','outerposition',[0 0 1 1]);
    image = imread(path);
    imagesc(image);
    O = ginput(2);
    result = inputdlg('Calibrate: ','Plot limits', 1,{'[ 0 0 1 1 ]'});
    result = eval(result{1});
    M = reshape( result, 2, 2)';
    %figure (1);

    ax = (M(2,1) - M(1,1)) / (O(2,1) - O(1,1));
    bx = M(1,1) - ax*O(1,1);
    ay = (M(2,2) - M(1,2)) / (O(2,2) - O(1,2));
    by = M(1,2) - ay*O(1,2);

    xy = [];
    series = {};
    n = 0;
    m = 0;
    hold on;

    p = plot(0,0);
    looping = 1;

    while looping
        [xi,yi,but] = ginput(1);

        if but == 1
            n = n+1;
            xy(:,n) = [xi;yi];
        elseif but == 3
            if n > 0
                xy(:,n) = [];
                n = n-1;
            end
        elseif but == 2
            % a key was pressed
            if uint8(get(fh,'CurrentCharacter')) == 13
                m = m+1;
                series{m} = xy';
                %disp(xy');
                xy = [];
                n=0;
            end  
        elseif isempty(but)
            looping = 0;
        end

        %clf;
        %hold off;
        %imagesc(image);
        %hold on;
        delete (p);
        if ~isempty(xy)
            p = plot(xy(1,:),xy(2,:),'ro');
        else
            p = plot(0,0);
        end
    end

    delete([path '.dat']);

    hold on;

    for i = 1:length(series)
        xy = series{i};
        plot (xy(:,1),xy(:,2),'ro');

        %Interpolate with a spline curve and finer spacing.
        t = 1:size(xy,1);
        ts = 1: 0.1: size(xy,1);
        xys = spline(t,xy',ts);

        % Plot the interpolated curve.
        plot(xys(1,:),xys(2,:),'b-');


        xy(:,1) = ax*xy(:,1) + bx;
        xy(:,2) = ay*xy(:,2) + by;

        save([path '.dat'],'xy','-ascii', '-tabs','-append');
        fid=fopen([path '.dat'],'a');
        fprintf(fid,'\n');
        fclose(fid);
    end
    %hold off
end