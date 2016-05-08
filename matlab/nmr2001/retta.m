%RETTA: funzione di una  retta. 
%USO:
%>> z = retta(a_b,const_data)
%INPUT: a_b: parametri della retta.
%  const_data: matrice di 1,2, o 3 righe di dati costanti.
% se no_righe==1, "z" sono le ordinate corrispondenti a "const_data" per
%						la trasformazione lineare; 
% se no_righe==3, "const_data"==[x_i;y_i;err_i] (ascisse, ordinate, ed 
%  errori sperimentali) e la funzione calcola il chi2;
% se no_righe==2, come sopra, con err_1=1;

function ff= retta(P,D);  	
theo= polyval(P,D(1,:)); 
if (size(D,1)==1), 
  ff=theo;
elseif (size(D,1)==2), 
  ff = sum((D(2,:) - theo).^2);
else
  ff = sum(((D(2,:) - theo)./D(3,:)).^2);
end


