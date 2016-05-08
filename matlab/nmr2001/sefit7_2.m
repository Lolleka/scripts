function f=sefit7_2(p,d);
%%Stim. spin echo per I=7/2, rilassamento Redfield. 
%%par=[T_long,Aoverall,A0 ]
nc=floor(length(p)/2);
p_=zeros(1,nc*8+1);
for k=0:nc-1;
  p_(k*8+(1:8)) = [p(k*2+1)./[1,6,15,28],...
		   p(k*2+2)*[1/84 3/44 150/728 1225/1716]];
end
p_(end)=p(end);
f=t2fitm(p_,d);


