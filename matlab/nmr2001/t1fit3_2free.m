function f=t1fit3_2free(p,d);
%%SR per I=3/2, rilassamento Redfield 'unconstrained'. 
%%par=[T_long,A_long,A_short,A0]
nc=floor(length(p)/3);
p_=zeros(1,nc*4+1);
for k=0:nc-1;
  p_(k*4+(1:4)) = [p(k*2+1)./[1,6], p(k*2+2:3)];
end
p_(end)=p(end);
f=t1fitm(p_,d);




