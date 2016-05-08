fn = 'Hsweep_32MHz.pna';
%fn = 'Fsweep_8T.pna';
%savedata = Integrate({fn, '-T'}); %, '-AMAGN'});
savedata = Integrate({fn, '-T', '-AMAGN'});
savedata = sortrows(savedata,1);
save(strrep(fn,'.pna','.dat'),'savedata','-ascii');
savedata = savedata(find(savedata(:,3) == 1),:);
savedata = savedata(find(savedata(:,2)>=0),:);
plot(savedata(:,1),savedata(:,2));