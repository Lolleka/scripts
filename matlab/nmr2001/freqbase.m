function f=freqbase(zn)
global NMRseq NMRshift
zn=round(zn(1));
if ((zn > size(NMRseq,2)) | (zn < 1))
  error('There''s no such zone')
end
f=NMRseq(3,zn) - 1e-6*(NMRshift(zn)+NMRseq(8,zn)*(...
    -NMRseq(7,zn)/2:NMRseq(7,zn)/2-1));


