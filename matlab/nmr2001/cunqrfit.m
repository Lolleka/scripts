function y=cunqrfit(p,dat);
n_c=floor(length(p)/3);
isoffs=(rem(length(p),3)>0);
gp=zeros(1,n_c*6+isoffs);
for ii=1:n_c;
  gp(ii*6+(-5:-3))=p(ii*3+(-2:0));
  gp(ii*6+(-2:-0))=p(ii*3+(-2:0)).*[.44718,.92417,.92417];
end;
y=gmfit(gp,dat);
