Hsweeps_freq = [16 20 23 25 32];
Hsweeps = {'Hsweep_16MHz.dat' 'Hsweep_20MHz.dat' 'Hsweep_23MHz.dat' 'Hsweep_25MHz.dat' 'Hsweep_32MHz.dat'};

figure(1);
p = []
for i = 1:length(Hsweeps)
    d = load(Hsweeps{i});
    bs = find(d(:,4)==0);
    d(bs,3) = 0;
    d = sortrows(d,1);
    d(:,3) = d(:,3)/max(d(:,3));
    a = 3*(max(d(:,3))-min(d(:,3)));
    p = [p plot(d(:,1),(d(:,3)-min(d(:,3)))*a+Hsweeps_freq(i))];
    hold on;   
    plot([0 10],[Hsweeps_freq(i) Hsweeps_freq(i)],'--','Color','green');
end
hold off;

Fsweeps_field = [1.8 5 5.7 8];
Fsweeps = {'Fsweep_1.8T.dat' 'Fsweep_5T.dat' 'Fsweep_5.7T.dat' 'Fsweep_8T.dat'};
%figure(2);
%a = 1;
%for i = 1:length(Fsweeps)
%    d = load(Fsweeps{i});
%    bs = find(d(:,3)==0);
%    d(bs,2) = 0;
%    d = sortrows(d,1);
%    d(:,3) = d(:,3)/max(d(:,3));
%    a = (max(d(:,3))-min(d(:,3)));
%    plot(d(:,1),(d(:,3)-min(d(:,3)))*a+Fsweeps_field(i));
%    hold on;
%    plot([0 40],[Fsweeps_field(i) Fsweeps_field(i)],'--','Color','green');
%end
%hold off;

figure(1);
hold on;
for i = 1:length(Fsweeps)
    d = load(Fsweeps{i});
    bs = find(d(:,4)==0);
    d(bs,3) = 0;
    d = sortrows(d,2);
    d(:,3) = d(:,3)/max(d(:,3));
    a = 0.7*(max(d(:,3))-min(d(:,3)));
    p = [p plot((d(:,3)-min(d(:,3)))*a+Fsweeps_field(i),d(:,2),'Color','red')];
    plot([Fsweeps_field(i) Fsweeps_field(i)], [10 35],'--','Color','green');
end
hold off;

uistack(p, 'top');