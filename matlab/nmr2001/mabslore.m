function [varargout]=mabslore(varargin);
%% wrapper of PEAKFIT: fit to many abs-lorentzians
pdescr=struct('fcn','almfit','np',3,'pnam',['Amp   ';'Cent  ';'Lambda'],...
	      'pty','acs');
varargout=cell(1,max(nargout,1));
[varargout{:}]=peakfit(pdescr,varargin{:});
