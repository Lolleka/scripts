%NUZQPLOT	Run-time plot of a quadrupole powder pattern within FMINUIT.
%		This function is in charge of drawing a run-time plot of the
%		fit status during the execution of fminuit, to which it is 
%		supplied as auxiliary function. It is called by ZEEMNUQ
%		and friends (znuqmgau, ...) whenever a 'CALL 5' minuit 
%		command is issued.
%
%Syntax:
%
%>> NUZQPLOT(PARS,DATA,FUNNAME,ERRMATRIX)
%
%INPUT ARGUMENTS
%PARS:    vector variational parameters;
%DATA:    matrix of constant data ([x_i;y_i], or [x_i;y_i;Dy_i], or 
%         [x_i;y_i;up_Dy_i;low_Dy]) arranged row-wise;  
%FUNNAME: name of a fit function [y,x]=funname(PARS,DATA), returning 2 row 
%	 vectors (the theoretical function y(x) and the independent variable
%	 x) if DATA is a single row.
%ERRMATRIX: unused.
%
%Never use this function directly. It is intended to be called only from within
%specialized fitting programs based on the POWDER package, like e.g. ZEEMNUQ, 
%ZNUQMGAU, ZNUQMBIGAU.
%
%
%   NMR2001: nuzqplot 
%   Revision: 1.0,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
%
%
%
%
%
%
function dummy=nuzqplot(par,aux,fcn,erma);
%%NULQPLOT: run-time plot during execution of gaussnuq (and friends).
global NMRplotsty

if isstruct(NMRplotsty),
  sym=NMRplotsty.sym;
else
  sym='*';
end

stat = 0;
if (nargin >2), 	
  if (isstr(fcn)),
    eval(sprintf('[y x]=%s(par,aux(1,:));',fcn));
    if (isvec(y) & isvec(x) & length(x)==length(y)),
      stat = 1;
    end
  end
end
if (stat),
  plot(aux(1,:),aux(2,:),sym,x,y,'b');
  drawnow;
else 
  disp('Nothing to plot');
end;
dummy=[];

%%%local macros
function r=isvec(a);
r=isnumeric(a) & (min(size(a))==1);
