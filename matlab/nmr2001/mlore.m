%MLORE		Fit of a spectrum to many normalized lorentzian peaks 
%		by means of FMINUIT, and plot. 
%		Each peak depends on three parameters: [area, center, width].
%		This program is a wrapper of PEAKFIT, to which it passes the 
%		appropriate function descriptor and all its input arguments. 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = MLORE(...)
%
%INPUT ARGUMENTS
%Input arguments (variable list) are optional; see PEAKFIT for details. 
%The default behaviour is one lorentzian peak, experimental data loaded from 
%the environment (see ENVIRONMENT).
%
%OUTPUT arguments are as-returned by FMINUIT, except CHISQUARE which is 
%normalized.
%
%
%EXAMPLES
%
%Fit of a spectrum with 2 lorentzian peaks (data loaded from the environment 
%and corrected for omega^2):
%
%>> [a b c] = mlore(2);
%
%or (an initial guess is specified)
%
%>> [a b c] = mlore([10 30 2 20 50 1.5]);
%
%Fit of a spectrum with 2 lorentzian peaks plus a constant offset (data loaded 
%from the environment and corrected for omega^2):
%
%>> [a b c] = mlore(2.5);
%
%or (an initial guess is specified)
%
%>> [a b c] = mlore([10 30 2 20 50 1.5 .03]);
%
%
%Fit of simulated data, inlined through the input list, with 2 lorentzians:
%
%>>x=0:100;
%>>y=10*exp(-.5*((x-30)/2).^2) + 20*exp(-.5*((x-30)/2).^2) + ...
%>>			      .05*(rand(size(x))-.5);
%>>dy=.02*ones(size(t));
%
%>> [a b c] = mlore(x,y,dy,2);
%or
%>> [a b c] = mlore([x;y;dy],2);
%
%or (an initial guess is specified)
%
%>> [a b c] = mlore(x,y,dy,[4 55 3 4 35 2]);
%or
%>> [a b c] = mlore([x;y;dy],[4 55 3 4 35 2]);
%
%
%   NMR2001: mlore
%   Revision: 1.0,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
%
%
function [varargout]=mlore(varargin);
%% wrapper of PEAKFIT: fit to many lorentzians
pdescr=struct('fcn','lmfit','np',3,'pnam',['Amp   ';'Cent  ';'Lambda'],...
	      'pty','acs');
varargout=cell(1,max(nargout,1));
[varargout{:}]=peakfit(pdescr,varargin{:});
