function [varargout]=mgaugau(varargin);
%% wrapper of PEAKFIT: fit to many concentrical bigaussians
pdescr=struct('fcn','mggfit','np',5,'pnam',['Aa   ';'Cent ';'Sa   ';...
		    'Ab/Aa';'Sb/Sa'],'pty','acsbb');
varargout=cell(1,max(nargout,1));
[varargout{:}]=peakfit(pdescr,varargin{:});
