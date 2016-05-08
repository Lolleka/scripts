function yax(a1,a2);
z=axis;
if (nargin>1), 
	if (isstr(a1)& ~isstr(a2)), z(4)=a2(1);
	elseif (isstr(a2) & ~isstr(a1)), z(3)=a1(1);
	else 	z(3:4)=[a1(1),a2(1)];
	end;
else
	z(3:4)=sort(a(1:2));
end
axis(z);
