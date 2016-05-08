function [x0, y0]=scrop(x_,y_);
%global NMRpar2d NMRamp;

if (length(x_)<2),
  x0=x_;
  y0=y_;
  return
end
[u,tinx]=sort(x_); 
mt_i=[x_(tinx);y_(tinx)].';
w=diff(mt_i(:,1)); 
inx=zeros(size(tinx));
cnt=1; bg=1; ed=1; 

for ii=1:length(tinx)-1;
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
if (~w(length(tinx)-1))  %multiplet
   [xx,yy]=max(mt_i(bg:ed,2)); 
   inx(cnt)=tinx(bg+yy-1);
else  %singlet
  inx(cnt)=tinx(bg);
end
  
%%keyboard
x0=x_(inx(1:cnt));
y0=y_(inx(1:cnt));
%zselect(inx(1:cnt));



