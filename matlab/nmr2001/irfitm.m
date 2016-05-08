function f=irfitm(p,d)
nc=(length(p)-1)/2;
f=ones(1,size(d,2))*p(end);

for k=1:nc;
    f=f+p(k+nc)*(1-2*exp(-(d(1,:)/p(k)))); 
end
if (size(d,1) == 2)
    f=sum((f-d(2,:)) .^2);
elseif (size(d,1) > 2)
    f=sum(((f-d(2,:))./d(3,:)).^2);
end
