function y=hygausswin(N,cent,wid,halfFL);
%%apodization windows for fid filtering
x=1:N;
y=exp(-.5*((x-cent)/wid).^4);
if (halfFL),
  y=y.*(x >= cent);
end
