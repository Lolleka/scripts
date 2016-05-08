%yflo: pone l'estremo inferiore dell'asse Y di un grafico =0
%USO: 
%>> yflo
function yflo(dummy);
z=axis; z(3)=0; axis(z); 
