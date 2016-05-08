function f=mbgfit(p,d);
%%multi bigaussian; each bigussian is made of 2 gauss with the same center

npk=floor(length(p)/5);
flg=(rem(length(p),5)>0); %%offset
pp=ones(1,6*npk+flg);
if (flg),
  pp(6*npk+1)=p(5*npk+1);
end;
for ii=1:npk;
  pp((1:6)+(ii-1)*6)=p([1:4,2,5]+(ii-1)*5);
end;
f=gmfit(pp,d);

