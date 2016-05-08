function [varargout]=srexp3_2(varargin);
%% wrapper of RELAXFIT: fit to magnetic Redfield SR, I=3/2
global NMRfitflags; 
fitflags_=NMRfitflags;
NMRfitflags.offs = 1;  %release ofset as a free fitting parameter
rdescr=struct('fcn','t1fit3_2','np',2,'offs',1,'interl',1,'pnam',...
	      ['T1 ';'Amp'],'pty', 'ta');
%%rdescr.fcn: fitting function
%%rdescr.np: No. of parameters for each component
%%rdescr.offs: %% %T: use offset
%%rdescr.interl: ==1 -> interlaced parameters, e.g.: T2#1, T2#2, A#1, A#2
%%rdescr.pnam: parameter names for each component
%%rdescr.pty: pameter type ('a'=amplitude, 'r'==rate, 't'=time constant, 
%%            'b'==beta

varargout=cell(1,max(nargout,1));
[varargout{:}]=relaxfit(rdescr,varargin{:});
NMRfitflags=fitflags_; %%restore 
