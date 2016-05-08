function [varargout]=seexp5_2(varargin);
%% wrapper of RELAXFIT: fit to magnetic Redfield SSE, I=5/2
rdescr=struct('fcn','sefit5_2','np',2,'offs',1,'interl',1,'pnam',...
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
