%MBIGAU		Fit of a spectrum to many concentric double gaussian peaks, 
%		and plot.
%		Each composite peak is made of 2 normalized gaussian curves 
%		sharing a common center, and depends on five parameters: 
%		[area#1, center, width#1, area#2, width#2]. 
%		This program is a wrapper of PEAKFIT: it builds the 
%		appropriate function descriptor and passes it to PEAKFIT, 
%		followed by all its input arguments. 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = MBIGAU(...)
%
%INPUT ARGUMENTS
%Input arguments (variable list) are optional, see PEAKFIT for details.
%The default behaviour is one bi-gaussian peak, experimental data loaded from 
%the environment (see ENVIRONMENT).
%
%OUTPUT arguments are as-returned by FMINUIT, except CHISQUARE which is 
%normalized.
%
%
%EXAMPLES
%
%Fit of a spectrum with 2 bi-gaussian peaks (data loaded from the environment 
%and corrected for omega^2):
%
%>> [a b c] = mbigau(2);
%
%or (an initial guess is specified)
%
%>> [a b c] = mbigau([10 30 1 15 3 20 80 1.5 25 4]);
%
%Fit of a spectrum with 2 bi-gaussian peaks plus a constant offset (data loaded 
%from the environment and corrected for omega^2):
%
%>> [a b c] = mbigau(2.5);
%
%or (an initial guess is specified)
%
%>> [a b c] = mbigau([10 30 1 15 3 20 80 1.5 25 4 .03]);
%
%
%Fit of simulated data, inlined through the input list, with 1 bi-gaussian:
%
%>>x=0:100;
%>>y=10*exp(-.5*((x-40)/2).^2) + 20*exp(-.5*((x-40)/5).^2) + ...
%>>			      .05*(rand(size(x))-.5);
%>>dy=.02*ones(size(t));
%
%>> [a b c] = mbigau(x,y,dy);
%or
%>> [a b c] = mbigau([x;y;dy]);
%
%or (an initial guess is specified)
%
%>> [a b c] = mbigau(x,y,dy,[10 55 3 20 8]);
%or
%>> [a b c] = mbigau([x;y;dy],[10 55 3 20 8]);
%
%
%   NMR2001: mbigau
%   Revision: 1.0,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
%
%
function [varargout]=mbigau(varargin);
%% wrapper of PEAKFIT: fit to many concentrical bigaussians
pdescr=struct('fcn','mbgfit','np',5,'pnam',['Amp_a  ';'Cent   ';'Sigma_a';...
		    'Amp_b  ';'Sigma_b'],'pty','acsas');
varargout=cell(1,max(nargout,1));
[varargout{:}]=peakfit(pdescr,varargin{:});
