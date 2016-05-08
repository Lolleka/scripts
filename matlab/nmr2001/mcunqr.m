%MCUNQR		Fit of Cu NQR spectra to many normalized gaussian 63Cu peaks, 
%		plus the corresponding 65Cu peak scaled in frequency and 
%		amplitude by the quadrupole moment ratio and the relative 
%		abundance of the two isotopes. 
%		Each peak component depends on three parameters: 
%		[area, center, width] of the 63Cu peak. 
%		This program is a wrapper of PEAKFIT: it builds the 
%		appropriate function descriptor and passes it to PEAKFIT, 
%		followed by all its input arguments. 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = MCUNQR(...)
%
%INPUT ARGUMENTS
%Input arguments (variable list) are optional, see PEAKFIT for details.
%The default behaviour is one gaussian 63Cu peak + its isotopic replica, 
%experimental data loaded from the environment (see ENVIRONMENT).
%
%OUTPUT arguments are as-returned by FMINUIT, except CHISQUARE which is 
%normalized.
%
%
%EXAMPLES
%
%Fit of a spectrum with two 63Cu gaussian peaks (data loaded from the 
%environment and corrected for omega^2):
%
%>> [a b c] = mcunqr(2);
%
%or (an initial guess is specified)
%
%>> [a b c] = mcunqr([10 30 2 20 50 1.5]);
%
%Fit of a spectrum with two 63Cu gaussian peaks plus a constant offset (data 
%loaded from the environment and corrected for omega^2):
%
%>> [a b c] = mcunqr(2.5);
%
%or (an initial guess is specified)
%
%>> [a b c] = mcunqr([10 30 2 20 50 1.5 .03]);
%
%
%   NMR2001: mcunqr
%   Revision: 1.0,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
%
%
function [varargout]=mcunqr(varargin);
%function [cf,errs,chi,emt]=mcunqr(varargin);
%% wrapper of PEAKFIT: fit to many 63Cu gaussian lines, + as many 65Cu replicas
pdescr=struct('fcn','cunqrfit','np',3,'pnam',['  Amp63_';' Cent63_';'Sigma63_'],...
  'pty','acs');
varargout=cell(1,max(nargout,1));

[varargout{:}]=peakfit(pdescr,varargin{:});
