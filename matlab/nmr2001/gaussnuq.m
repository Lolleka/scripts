%GAUSSNUQ	Fitting and plot of spectral data to a powder quadrupole 
%		pattern with gaussian broadening for B>>Q. 
%		This function is obsolete: call ZEEMNUQ istead.
%
%   NMR2001: gaussnuq 
%   Revision: 1.1,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
function [varargout]=gaussnuq(varargin);

warning(['GAUSSNUQ is obsolete and may be removed in the future. Use ' ...
	 'ZEEMNUQ instead']);
varargout=cell(1,max(nargout,1));
[varargout{:}]=zeemnuq(varargin{:});

%% Copyright (C) G. Allodi
%% Revision 2.0   14-sept-2002
