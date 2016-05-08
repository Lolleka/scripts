%erbawid: sets errorbar's width to a user defined value WIDTH
%USAGE: erbawid(handle,WIDTH);
%G.A. scripsit

function y=erbawid(id,wid);

for ii=1:length(id); 
   if strcmp(get(id(ii),'type'),'line'); %if handle to lines...
	x=get(id(ii),'xdata');
	n=length(x);
	if isnan(x(-6 +(9:9:n))), %true for errobars
		x(-5 +(9:9:n))=x(-8 +(9:9:n))-wid;
		x(-2 +(9:9:n))=x(-8 +(9:9:n))-wid;
		x(-4 +(9:9:n))=x(-8 +(9:9:n))+wid;
		x(-1 +(9:9:n))=x(-8 +(9:9:n))+wid;
	end
	set(id(ii),'xdata',x);
   end
end
