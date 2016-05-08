function mt=readdata(fnam,vmod);

if (nargin<2), %vmod, print header
	vmod=0;
else 
	vmod=max(abs(vmod));
end

maxrow=0; maxclmn=0; rowcnt=0;

id=fopen(fnam,'r');
if (id <= 0),
  error(sprintf('File %s could not be opened',fnam));
end;

while 1;
  if (feof(id)), %check after EOF. NB: feof is non 0 at EOF even b4 attempting 
	break;   %to read. This differs from C's feof behavior
  end;
  str=fgets(id); %read str: 
  rowcnt=rowcnt+1;
  [vec,n]=sscanf(str,'%g ');
  if (n>1),
	maxrow=maxrow+1;
        maxclmn=max(maxclmn,n);	
  end;
end;

if (fseek(id,0,'bof'))
	str=sprintf('File %s could not be rewound',fnam);  
	error(str);
end;

mt=zeros(maxrow,maxclmn);
cnt_row_ok=0;

for ii=1:rowcnt;
  if (feof(id)),
	error(sprintf('Unexpected EOF in %s',fnam));
  end;
  str=fgets(id); %read str
  [vec,n]=sscanf(str,'%g ');

  if (n>1),
	cnt_row_ok=cnt_row_ok+1;
	mt(cnt_row_ok,1:n)=vec.';
%%%
	mt(cnt_row_ok,n+1:maxclmn) = (n+1:maxclmn)*NaN;
					%a shorter row is padded with NANs
%%%
	vmod=0; 
  else 
	if (vmod),
		disp(str(1:length(str)-1));
	end;	
  end;
end;

mt=mt(1:cnt_row_ok,:);

fclose(id);

