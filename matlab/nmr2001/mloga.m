%MLOGA		Fit of a spectrum to many gaussian-or-lorentzian peaks, 
%		and plot. 
%		Each peak depends on three parameters: [area, center, width].
%		A peak is gaussian if its width is positive, lorentzian 
%		if width is negative. In the latter case the peak amplitude 
%		and width are -area and -width, respectively.
%		This program is a wrapper of PEAKFIT, to which it passes the 
%		appropriate function descriptor and all its input arguments. 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = MLOGA(...)
%
%INPUT ARGUMENTS
%Input arguments (variable list) are optional, see PEAKFIT for details.
%The default behaviour is one gaussian peak, experimental data loaded from 
%the environment (see ENVIRONMENT).
%
%OUTPUT arguments are as-returned by FMINUIT, except CHISQUARE which is 
%normalized.
%
%NB: since the peak shape depends on it parameter, the fit may switch from
%    gaussian to lorentzian peaks, and vice versa.
%
%
%EXAMPLES
%
%Fit of a spectrum with 2 gaussian peaks (data loaded from the environment 
%and corrected for omega^2):
%
%>> [a b c] = mloga(2);
%
%Fit of a spectrum with 1 gaussian + 1 lorenzian peak + a constant offset
%(NB: the peak shapes are kept by the fit if FMINUIT returns a positive and 
%a negative width, respectively). 
%Data loaded from the environment and corrected for omega^2; an initial guess 
%is specified.  
%
%>> [a b c] = mloga([10 30 2 -20 50 -1.5 .06]);
%
%
%Fit of simulated data, inlined through the input list, with 1 gaussian + 
%1 lorenzian peak; an initial guess is specified:
%
%>>x=0:100;
%>>y=40*exp(-.5*((x-30)/2).^2) +60*3.0./((x-55).^2 + 3^2)/pi + ...
%>>			      .05*(rand(size(x))-.5);
%>>dy=.02*ones(size(t));
%
%>> [a b c] = mloga(x,y,dy,[44 38 1 -77 69 -4]);
%or
%>> [a b c] = mloga([x;y;dy],[44 38 1 -77 69 -4]);
%
%
%   NMR2001: mloga
%   Revision: 1.0,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
%
%
function [varargout]=mloga(varargin);
%% Wrapper of PEAKFIT: fit to many gaussian or lorentian peaks.
%% The lineshape is swiched to lorentzian by reversing together the sign of 
%%peak amplitude and width.
pdescr=struct('fcn','lgmfit','np',3,'pnam',['Amp  ';'Cent ';'Sigma'],...
	      'pty','acs');
varargout=cell(1,max(nargout,1));

[varargout{:}]=peakfit(pdescr,varargin{:});
