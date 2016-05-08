%merges 1 spectrum by keeping the max value for each set of "equal" experiments
%usage [mt_o, inx]=selfmerge(mt_i);
% mt_* =[parameter(:),values(:),...]

function [mt_o, inx]=selfmerge(mt_i,mod);

if (nargin < 2), mod=1; end;

[u,tinx]=sort(mt_i(:,1)); mt_i = mt_i(tinx,:);
w=diff(mt_i(:,1)); 
mt_o=zeros(size(mt_i)); 
inx=zeros(size(tinx));

cnt=0; bg=1; ed=1; 

for ii=1:length(tinx)-1;
	if (~w(ii)),
		ed=ed+1; 
	else
		cnt=cnt+1; 
		[xx,yy]=max(mt_i(bg:ed,2)); 
		if (~mod), 
			mt_o(cnt,:)=mean(mt_i(bg:ed,:)); 
		else
			mt_o(cnt,:)=mt_i(bg+yy-1,:); 
		end;
		inx(cnt)=tinx(bg+yy-1);
		bg=ed+1; ed=bg; 
	end;
end;

%store buffer remaining after for loop 

cnt=cnt+1; 
[xx,yy]=max(mt_i(bg:ed,2)); 
if (~mod), 
	mt_o(cnt,:)=mean(mt_i(bg:ed,:)); 
else
	mt_o(cnt,:)=mt_i(bg+yy-1,:); 
end;
inx(cnt)=tinx(bg+yy-1);

%trim final matrix
mt_o=mt_o(1:cnt,:); inx=inx(1:cnt);
