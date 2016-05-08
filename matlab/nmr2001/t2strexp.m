function [varargout]=t2strexp(varargin);
%% wrapper of RELAXFIT: fit to (many) stretched expo decay 
rdescr=struct('fcn','t2stfitm','np',3,'offs',1,'interl',1,'pnam',...
	      ['T2/2';'beta';'Amp '],'pty', 'tba');
%%rdescr.fcn: fitting function
%%rdescr.np: No. of parameters for each component
%%rdescr.offs: %% ==1 -> enable offset
%%rdescr.interl: ==1 -> interlaced parameters, e.g.: T2#1, T2#2, A#1, A#2
%%rdescr.pnam: parameter names for each component
%%rdescr.pty: pameter type ('a'=amplitude, 'r'==rate, 't'=time constant, 
%%            'b'==beta

varargout=cell(1,max(nargout,1));
[varargout{:}]=relaxfit(rdescr,varargin{:});
