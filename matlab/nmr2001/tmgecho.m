%TMGECHO  Analysis of TecMag NMR data: apodization and Fourier transform. 
% Usage: 
%  TMGECHO(FILENAME,...)
%
%(dots stand for a variable list of OPTION arguments). 
%TMGECHO is a wrapper of ECHOES; it is the same as  
%>> ECHOES(FILENAME,'-M0',...) 
%
%
%   NMR2001: tmgecho 
%   Revision: 2.0,  14-Sept-2002
%   Copyright (c) 2001-03 by Giuseppe Allodi
%
%
function tmgecho(varargin);
varargin{end+1}='-M0';
echoes(varargin{:});

