%T2MEXP		Fit and plot of relaxation data to a many-exponential decay 
%		+ constant offset, as for spin-spin relaxation (T2).
%		Fitting is performed by the FMINUIT minimization engine. 
%		The number N of exponential components is determined by the 
%		length L of the fitting parameter vector, L = 2*N+1. 
%		The fitting parameters are interlaced, i.e. they are oredered
%		like [Tau_1, ..., Tau_N, Amp_1,..., Amp_N, Amp_offset],
%		where Tau_k are time constants, and Amp_k are their respective 
%		amplitudes. The last parameter is the amplitude of a constant 
%		offset.	
%		This program is a wrapper of RELAXFIT, to which it passes the 
%		appropriate function descriptor and all its input arguments. 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = T2MEXP(...)
%
%INPUT ARGUMENTS
%Input arguments (variable list) are optional, see RELAXFIT for details.
%The default behaviour is one decay component, experimental data loaded from 
%the environment (see ENVIRONMENT).
%
%OUTPUT arguments are as-returned by FMINUIT, except CHISQUARE which is 
%normalized.
%
%
%EXAMPLES
%
%Fit of a T2 relaxation with 2 exponential components plus a constant offset 
%(data loaded from the environment):
%
%>> [a b c] = t2mexp(2);
%
%or (an initial guess is specified)
%
%>> [a b c] = t2mexp([1 .1 30 20 0]);
%
%Fit of simulated data, inlined through the input list, with 2 exponentials + 
%offset:
%
%>>x=0:100;
%>>y=10*exp(-x/18) + 20*exp(-x/38) + .05*(rand(size(x))-.5);
%>>dy=.02*ones(size(x));
%
%>> [a b c] = t2mexp(x,y,dy,2);
%or
%>> [a b c] = t2mexp([x;y;dy],2);
%
%or (an initial guess is specified)
%
%>> [a b c] = t2mexp(x,y,dy,[4 55 3 4 35 2]);
%or
%>> [a b c] = t2mexp([x;y;dy],[4 55 3 4 35 2]);
%
%   NMR2001: t2mexp 
%   Revision: 1.0,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
function [varargout]=t2mexp(varargin);
%% wrapper of RELAXFIT: fit to many expo-decay 
rdescr=struct('fcn','t2fitm','np',2,'offs',1,'interl',1,'pnam',...
	      ['T2/2';'Amp '],'pty', 'ta');
%%rdescr.fcn: fitting function
%%rdescr.np: No. of parameters for each component
%%rdescr.offs: %% %T: use offset
%%rdescr.interl: ==1 -> interlaced parameters, e.g.: T2#1, T2#2, A#1, A#2
%%rdescr.pnam: parameter names for each component
%%rdescr.pty: pameter type ('a'=amplitude, 'r'==rate, 't'=time constant, 
%%            'b'==beta

varargout=cell(1,max(nargout,1));
[varargout{:}]=relaxfit(rdescr,varargin{:});
