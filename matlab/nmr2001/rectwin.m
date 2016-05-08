function y=rectwin(N,cent,wid,halfFL);
%%apodization windows for fid filtering
x=1:N;
if (halfFL),
  y=((x-cent) <= wid) & ((x-cent) >= 0);
else
  y=(abs(x-cent) <= wid);
end

