function f=t1fit3_2(p,d)
nc=(length(p)-1)/2;
ep=zeros(1,4*nc+1);

for k=0:nc-1;
  ep(k*4+1:2)=p(k*2+1)*[1/6 1];
  ep(k*4+3:4)=p(k*2+2)*[.6 .4];
end
ep(nc*4+1)=p(2*nc+1);
f=t1fitm(ep,d);

  
