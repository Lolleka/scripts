%xzoom: espande l'asse X di un grafico per un un fattore x_fac (default =.1)
%USO: 
%>> xzoom(x_fac)

function f=xzoom(x);
if ~nargin, x=.1; end;
z=axis; z(2)=z(1)+(z(2)-z(1))*x(1); axis(z);
