function [varargout]=mabsgau(varargin);
%% wrapper of PEAKFIT: fit to many abs-gaussians
pdescr=struct('fcn','agmfit','np',3,'pnam',['Amp  ';'Cent ';'Sigma'],...
	      'pty','acs');
varargout=cell(1,max(nargout,1));
[varargout{:}]=peakfit(pdescr,varargin{:});
