function cl=cellify(str);
%convert a string containing '\n' into a cell array of strings 
if (~isstr(str)),
  cl={};
  return
end
ix=find(str==10);
n=length(ix);

if (~n),
  cl={str};
else
  cl=cell(n+1,1);
  cl{1}=str(1:ix(1)-1);
  cl{end}=str(ix(end)+1:end);
  for k=2:n;
    cl{k}=str(ix(k-1)+1:ix(k)-1);
  end
end
