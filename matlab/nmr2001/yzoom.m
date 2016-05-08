%yzoom: espande l'asse Y di un grafico per un un fattore y_fac (default =.1)
%USO: 
%>> yzoom(y_fac)

function f=yzoom(x);
if ~nargin, x=.1; end;
z=axis; z(4)=z(3)+(z(4)-z(3))*x(1); axis(z);
