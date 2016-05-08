function f=mggfit(p,d);
%%Multi bigaussian; each bigussian is made of 2 gauss with the same center.
%%Amplitude and width of 2nd component are experessed relative to those 
%%of the 1st (i.e., as ratio A2/A1, s2/s1)

npk=floor(length(p)/5);
flg=(rem(length(p),5)>0); %%offset
pp=ones(1,6*npk+flg);
if (flg),
  pp(6*npk+1)=p(5*npk+1);
end;
for ii=1:npk;
  v=p([1:5]+(ii-1)*5);
  pp((1:6)+(ii-1)*6) = [v(1:3) v(1)*v(4) v(2) v(3)*v(5)];
end;
f=gmfit(pp,d);

