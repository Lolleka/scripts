%SLRECHO  Analysis of Stelar NMR data: apodization and Fourier transform. 
% Usage: 
%  SLRECHO(FILENAME,...)
%
%(dots stand for a variable list of OPTION arguments). 
%SLRECHO is a wrapper of ECHOES; it is the same as  
%>> ECHOES(FILENAME,'-M1',...) 
%
%
%   NMR2001: slrecho 
%   Revision: 2.0,  14-Sept-2002
%   Copyright (c) 2001-03 by Giuseppe Allodi
%
%
function slrecho(varargin);
varargin{end+1}='-M1';
echoes(varargin{:});
