%MNPLOT		Run-time plot of a fit with FMINUIT.
%		This function is in charge of drawing a run-time plot of the
%		fit status during the execution of fminuit, to which it is 
%		supplied as auxiliary function. It is called by PEAKFIT, 
%		RELAXFIT, and their wrappers whenever a 'CALL 5' minuit 
%		command is issued.
%
%Usage:
%
%>> MNPLOT(PARS,DATA,FUNNAME,ERRMATRIX)
%
%INPUT ARGUMENTS
%PARS:    vector variational parameters;
%DATA:    matrix of constant data ([x_i;y_i], or [x_i;y_i;Dy_i], or 
%         [x_i;y_i;up_Dy_i;low_Dy]) arranged row-wise;  
%FUNNAME: name of function funname(PARS,DATA), so that it returns a row vector
%	 (the theoretical function values) if DATA is a single row.
%ERRMATRIX: unused.
%
%
%EXAMPLES
%
%>> fminuit('t2fitm','mnplot',[20 5 0],[0:100;5*exp(-(0:100)/20)]);
%  CALL 5
%
%produces the following Matlab call:
%
%>> mnplot([20 5 0],[0:100;5*exp(-(0:100)/20)],'t2fitm')
%
%
%   NMR2001: mnplot 
%   Revision: 1.0,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
%
%
%
%
%
%
function dummy=mnplot(par,aux,fcn,erma);
stat = 0;
if (nargin >2), 	
  if (isstr(fcn)),
    l=min(aux(1,:)); 
    r=max(aux(1,:)); 
    n=max(500,size(aux,2));
    x =l:(r-l)/n:r;
    y= feval(fcn,par,x);
    stat = (length(y) == length(x));
  end;
end;
if (stat),
  plot(aux(1,:),aux(2,:),'*',x,y,'-');
  drawnow;
else 
  disp('Nothing to plot');
end;
%if (nargin == 4), disp(erma); end;
dummy=[];
