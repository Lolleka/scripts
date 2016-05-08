function y=expowin(N,cent,wid,halfFL);
%%apodization windows for fid filtering
x=1:N;
y=exp(-abs((x-cent)/wid));
if (halfFL),
  y=y.*(x >= cent);
end
