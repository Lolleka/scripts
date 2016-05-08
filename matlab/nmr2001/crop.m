function crop;
global NMRpar2d NMRamp;


N=length(NMRpar2d.val);
if (N<2),
  return
end
[u,tinx]=sort(NMRpar2d.val); 
mt_i=[NMRpar2d.val(tinx);NMRamp(tinx)].';
w=diff(mt_i(:,1)); 
inx=zeros(1,N);
cnt=1; bg=1; ed=1; 

for ii=1:N-1;
  if (~w(ii)),
    ed=ed+1; 
  else
    [xx,yy]=max(mt_i(bg:ed,2)); 
    inx(cnt)=tinx(bg+yy-1);
    bg=ed+1; 
    ed=bg; 
    cnt=cnt+1; 
  end;
end

%termination
if (~w(N-1))  %multiplet
   [xx,yy]=max(mt_i(bg:ed,2)); 
   inx(cnt)=tinx(bg+yy-1);
else  %singlet
  inx(cnt)=tinx(bg);
end

zselect(inx(1:cnt)); %does everything


