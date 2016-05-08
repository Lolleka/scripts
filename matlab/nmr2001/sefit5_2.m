function f=sefit5_2(p,d);
%%Stim. spin echo per I=5/2, rilassamento Redfield. 
%%par=[T_long,A_overall,A0]
nc=floor(length(p)/2);
p_=zeros(1,nc*6+1);
for k=0:nc-1;
  p_(k*6+(1:6)) = [p(k*2+1)./[1,6,15], p(k*2+2)*[1/35 8/45 50/63]];
end
p_(end)=p(end);
f=t2fitm(p_,d);




