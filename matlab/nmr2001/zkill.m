%ZKILL: delete a number of 2-D zones (experiments) from workspace. 
%USAGE:
%  zkill(indices)
% kills zones of order indices ("indices" is a vector);
%  zkill(indices,1)
% kills the complement of indices

function dummy=zkill(ix,fl);
global TDfids FDfids NMRamp NMRseq NMRnoise NMRpar2d;
global NMRmisc;
matnam0={'TDfids','NMRseq'};
matnam1={'FDfids','NMRamp','NMRnoise','NMRpar2d.val'};

if (nargin <2), 
  fl=0; 
end
dim2d=size(NMRseq,2);
%validate indices
ix(find(ix < 1 | ix > dim2d))=[];
if (~length(ix)),
  disp('No valid index, nothing done');
  return
end
%check integrity of workspace
if (all(NMRseq(7,:)))
  matnam=[matnam0 matnam1];
elseif (any(NMRseq(7,:)))
  error('Analysis unsynced to time-domain data, cannot kill/select zones');
else
  matnam=matnam0;
end
  
for  ii=1:length(matnam);
  if (eval(['size(',matnam{ii},',2) ~= dim2d & length(',matnam{ii},') ~= 1 ' ...
        '| ~isnumeric(',matnam{ii},')'])),
    error('Corrupted workspace'); 
  end   
end

%do zone deletion
for  ii=1:length(matnam);
  if (eval(['size(',matnam{ii},',2) == dim2d'])),	  
    if (fl),
      eval([matnam{ii},'=',matnam{ii},'(:,ix);']);
    else	
      eval([matnam{ii},'(:,ix)=[];']);
    end	
  end
end

%NMRmisc has different semantics with different datasets
if (length(NMRmisc)),
  if (iscell(NMRmisc)) %ac300, .tnt 
      if (iscell(NMRmisc{1})) %ac300
          for k=1:length(NMRmisc);
              if (fl), %select
                  NMRmisc{k}=NMRmisc{k}(ix);
              else   %kill
                  NMRmisc{k}(ix) = [];
              end
          end
      else %.tnt
          if (fl), %select
              NMRmisc=NMRmisc(ix);
          else   %kill
              NMRmisc(ix) = [];
          end
      end
  else	%.pna, .sdf
    if (fl), %select
      NMRmisc=NMRmisc(:,ix);
    else   %kill
      NMRmisc(:,ix) = [];
    end
  end	
end
