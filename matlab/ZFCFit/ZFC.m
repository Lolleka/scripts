files = dir ('.\Tutti50\*.dat');

figure(1);

TCut1 = [30 115 50];
TCut2 = [160 280 200];
for i=2:2%1:length(files)
    files(i).name = ['.\Tutti50\' files(i).name];
    [x res] = ZFCFit(files(i).name,TCut1(i),TCut2(i),1E6,2,1E-30,[90 900 800]);
    
    data = load(files(i).name);
    plot (data(:,1),data(:,2));
    hold on;
    
    plot(res(1,:),res(2,:),'Color','red');
    hold on;
    
    savedata = res';
    save([files(i).name '.fit'], 'savedata', '-ascii');
    savedata = x;
    save([files(i).name '.fitpar'], 'savedata', '-ascii');
    
    refresh(gcf);
end
hold off;