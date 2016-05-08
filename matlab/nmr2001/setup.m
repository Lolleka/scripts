global FDfids TDfids NMRseq NMRpath;
global NMRamp NMRshift NMRnoise NMRmisc;
global commands parnames stepbounds;
global NMRpar2d NMRplotsty NMRapod NMRfftwin NMRfitflags
global NMRstat


%defaults
NMRpath = '';
%defaults - initialize cell/struct variables
NMRapod = struct('typ','hygausswin','fwid',16e4,'half',0);
NMRfftwin = struct('foffs',0,'fwid',0,'funit',1); 
                                         %at center of FFT, use channel units
NMRfitflags = struct('batch',0,'offs',0);
NMRplotsty = struct('sym','o','line','-m','ebar',1,'drawpeaks',0,'grid',...
		    0,'ylog',0);
NMRpar2d = struct('name2d','none','val',[]);
NMRstat={}; 

if (exist('config.mat')),
  load config.mat
end
format short e







