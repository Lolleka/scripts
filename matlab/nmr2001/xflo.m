%xflo: pone l'estremo inferiore dell'asse X di un grafico =0
%USO: 
%>> xflo
function xflo(dummy);
z=axis; z(1)=0; axis(z); 
