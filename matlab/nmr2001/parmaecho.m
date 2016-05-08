%PARMAECHO  Analysis of pna data: apodization and Fourier transform. 
% Usage: 
%  PARMAECHO(FILENAME,...)
%
%(dots stand for a variable list of OPTION arguments). 
%PARMAECHO is a wrapper of ECHOES; it is the same as  
%>> ECHOES(FILENAME,'-M2',...) 
%
%
%   NMR2001: parmaecho 
%   Revision: 1.1,  18-Apr-2005
%   Copyright (c) 2001-05 by Giuseppe Allodi
%
%
function parmaecho(varargin);
varargin{end+1}='-M2';
echoes(varargin{:});
