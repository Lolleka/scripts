%TELL2D (private): 
%determines the type of experiment and returns the relevant 2D par
function [par2d,name2d,varseqFL]=tell2d(fd); %fd=1: freq. domain 
global NMRseq NMRmisc;

if (~nargin), %default: time domain
  fd=0;
end

par2d=NMRseq(3,:); %init value=frequency 1st zone
name2d='Frequency';
found=0;
if (size(NMRseq,2) > 1), 
  for ii=[3,9:size(NMRseq,1)];
    if any(diff(NMRseq(ii,:))),
      if (ii>4),
	par2d=NMRseq(ii,:);
	name2d=sprintf('Delay #%d',ii-8);
      end
      found=1;
      break
    end
  end
  % on Stelar and Hyrespect, seek NMRmisc as well
  if (~found & length(NMRmisc)),
    if (isnumeric(NMRmisc)),
      misc=NMRmisc;
    elseif (iscell(NMRmisc) & isnumeric(NMRmisc{end})),
      misc=NMRmisc{end}; 
    else 
      misc=[];
    end
    if (length(misc)),
      for ii=1:size(misc,1);
        if any(diff(misc(ii,:))),
          par2d=misc(ii,:);
          name2d=sprintf('Misc #%d',ii);
          break;
        end
      end
    end
  end
  if (fd),
    rng=7:8;
  else
    rng=[1:2 5];
  end
  varseqFL=any(diff(NMRseq(rng,:)')); 
  %True if acq parameters differ throughout 2D experiment
else
  varseqFL = 0;
end

