Hsweeps_freq = [16 20 23 25 32];
Hsweeps = {'Hsweep_16MHz.dat' 'Hsweep_20MHz_bis.dat' 'Hsweep_23MHz.dat' 'Hsweep_25MHz.dat' 'Hsweep_32MHz.dat'};

figure(1);
a = 3;
p = []
for i = 1:length(Hsweeps)
    d = load(Hsweeps{i});
    bs = find(d(:,3)==0);
    d(bs,2) = 0;
    d = sortrows(d,1);
    d(:,2) = d(:,2)/max(d(:,2));
    a = 2*(max(d(:,2))-min(d(:,2)));
    p = [p plot(d(:,1),(d(:,2)-min(d(:,2)))*a+Hsweeps_freq(i))];
    hold on;
    plot([0 10],[Hsweeps_freq(i) Hsweeps_freq(i)],'--','Color','green');
end
hold off;

Fsweeps_field = [1.8 5 5.7 8];
Fsweeps = {'Fsweep_1.8T.dat' 'Fsweep_5T.dat' 'Fsweep_5.7T.dat' 'Fsweep_8T.dat'};

figure(2);
a = 1;
for i = 1:length(Fsweeps)
    d = load(Fsweeps{i});
    bs = find(d(:,3)==0);
    d(bs,2) = 0;
    d = sortrows(d,1);
    d(:,2) = d(:,2)/max(d(:,2));
    a = (max(d(:,2))-min(d(:,2)));
    plot(d(:,1),(d(:,2)-min(d(:,2)))*a+Fsweeps_field(i));
    hold on;
    plot([0 40],[Fsweeps_field(i) Fsweeps_field(i)],'--','Color','green');
end
hold off;



figure(1);
hold on;
a = 0.5;
for i = 1:length(Fsweeps)
    d = load(Fsweeps{i});
    bs = find(d(:,3)==0);
    d(bs,2) = 0;
    d = sortrows(d,1);
    d(:,2) = d(:,2)/max(d(:,2));
    a = (max(d(:,2))-min(d(:,2)));
    p = [p plot((d(:,2)-min(d(:,2)))*a+Fsweeps_field(i),d(:,1),'Color','red')];
    plot([Fsweeps_field(i) Fsweeps_field(i)], [10 35],'--','Color','green');
end
hold off;

uistack(p, 'top');