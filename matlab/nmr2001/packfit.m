function  f=packfit(mat,a,b,c);
%%PACKFIT: appende il risultato di un fit (da fminuit) ad una matrice, in 
%%formato impaccato [...; dato#1,errore#1,...,dato#n,errore#n,chi2]
%%USO: 
%%>> out_matrix = packfit(in_matrix, NMRdata1, errors,chi2]
if (nargin < 3), 
  error('3 args at least');
elseif (nargin==3), 
  c=b; b=a; a=mat; mat=[];
end

rw= zeros(size(a,1),2*size(a,2)+1);
rw(:,1:2:2*size(a,2))=a;
rw(:,2:2:2*size(a,2))=b;
rw(:,end)=c;
f=[mat;rw];


%function f=isvector(v)
%f=min(size(v))==1;
