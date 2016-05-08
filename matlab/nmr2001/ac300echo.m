%AC300ECHO  Analysis of Bruker AC300 data: apodization and Fourier transform. 
% Usage: 
%  AC300ECHO(FILENAME,...)
%
%(dots stand for a variable list of OPTION arguments). 
%AC300ECHO is a wrapper of ECHOES; it is the same as  
%>> ECHOES(FILENAME,'-M0',...) 
%
%
%   NMR2001: ac300echo 
%   Revision: 1.0,  20-Apr-2010
%   Copyright (c) 2001-10 by Giuseppe Allodi
%
%
function ac300echo(varargin);
varargin{end+1}='-M3';
echoes(varargin{:});

