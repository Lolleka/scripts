function f=agmfit(p,d);
%Multi abs-complex-gaussian: theoretical function & chi^2
%Abs-complex-gaussian is defined as abs(g + i*HT(g)), where g is a
%gaussian (absorption) and HT(g) its Hilbert transform (dispersion)



npk=floor(length(p)/3);
flg=(rem(length(p),3)>0); %%offset
if (flg),
  f=p(3*npk+1);
else 
  f=0;
end;
for k=1:npk;
  p0=p((1:3)+(k-1)*3);
  x=d(1,:)-p0(2);
  f=f + abs(gmfit(p0,d(1,:))+...
            i*p0(1)*x.*kummerm(1,1.5,-.5*(x./p0(3)).^2)/(pi*p0(3)^2));
end;
%%f=abs(f);
if (size(d,1)>2)
  f=sum(((f-d(2,:))./d(3,:)).^2);
elseif(size(d,1)==2)
  f=sum((f-d(2,:)).^2);
end 



