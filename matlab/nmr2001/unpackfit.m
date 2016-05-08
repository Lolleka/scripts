function  [a,b,c]=unpackfit(mat);
%%UNPACK: spacchetta una matrice di dati di fitting in formato impaccato 
%%(cfr. PACKER) in tre matrici di: dati, errori, chi2
%%USO: 
%%>> [NMRdata1, errors,chi2] = unpackfit(packed_NMRdata1);

sz=size(mat);
a=zeros(sz(1),round(sz(2)/2-1)); 
b=a;
for ii = 1:(round(sz(2)/2-1)); 
  a(:,ii)=mat(:,ii*2-1); 
  b(:,ii)=mat(:,ii*2); 
end
c=mat(:,sz(2));


