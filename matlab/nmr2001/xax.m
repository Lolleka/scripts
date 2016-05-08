function xax(a1,a2);
z=axis;
if (nargin>1), 
	if (isstr(a1)& ~isstr(a2)), z(2)=a2(1);
	elseif (isstr(a2) & ~isstr(a1)), z(1)=a1(1);
	else 	z(1:2)=[a1(1),a2(1)];
	end;
else
	z(1:2)=sort(a(1:2));
end
axis(z);
