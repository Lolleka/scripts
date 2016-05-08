%RELAXFIT	Gateway program for fitting and plot of relaxation data.
%		Data are fitted by the FMINUIT minimization engine to a multi-
%		component model, specified by a suitable descriptor.  
%		NB.: This routine is not intended to be called directly by
%		a user. Call instead the appropriate wrapper function (e.g. 
%		t2mexp) which takes care of passing the descriptor of the 
%		appropriate fitting function.
% 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = RELAXFIT([DESCRIPTOR],...)
%
%(dots stand for a variable argument list). Optionally, the argument list
%may contain experimental DATA (by default, they are loaded from the 
%environment, see also ENVIRONMENT), initial GUESS of fitting parameters, 
%number of COMPONENTS, and option SWITCHES.
%
%
%INPUT PARAMETERS
%
%DESCRIPTOR is a struct specifying the fit function and its properties. The
%default function is a multi-exponential decay.
%
% DESCRIPTOR.fcn (string) 
%  The fit function name. The fit function has the prototype y=fcn(P,D), 
%  where P is a vector of fit parameters, and D is 1-,2-, or 3-row matrix 
%  of experimental data. 
%  The length of P is variable: 
%    length(P) = DESCRIPTOR.np * NC + DESCRIPTOR.offs
%  where NC is the number of components.
%  If D is a row-vector, it is interpreted as the time-base vector and fcn() 
%  returns the corresponding values of the theoretical function; otherwise
%  it is interpreted as [t;Y_t] or [t;y_t;Delta_Y_t] and the function returns 
%  the unnormalized Chi^2
%
% DESCRIPTOR.np (scalar) 
%  The number of fit parameters for each component. 
%  For instance, DESCRIPTOR.np=2 (amplitude, time constant) for a 
%  multiexponential model (T2FITM).  
%
% DESCRIPTOR.offs (boolean)
%  Flag: true if a constant offset is allowed for. In that case, the offset
%  is the last element of the fit parameter vector P. 
%
% DESCRIPTOR.interl (boolean)
%  Parameter ordering. DESCRIPTOR.interl is to be set true if fit parameters 
%  belonging to different components are interlaced. 
%  Let's assume, for instance, a 3-component model where each component 
%  depends on 2 parameters A_, B_: parameters are interlaced if they are 
%  ordered like
%    [A_1, A_2, A_3, B_1, B_2, B_3];
%  non-interlaced if they are ordered like
%    [A_1, B_1, A_2, B_2, A_3, B_3]; 
%
% DESCRIPTOR.pnam (matrix of characters)
%  Parameter names for a component. The number of rows of DESCRIPTOR.pnam 
%  equals DESCRIPTOR.np. Each row is a parameter name.
%
% DESCRIPTOR.pty (string) 
%  Parameter types. The length of DESCRIPTOR.pty equals DESCRIPTOR.np.
%  For each parameter of a component, a character specifies the parameter 
%  type. Characters are:
%  'a': amplitude
%  't': time constant;
%  'r': rate;
%  'b': stretching parameter Beta.
%
% 
%DATA 
%Data may be inlined as a 2-row [t_i;a_i], a 3-row [t_i;a_i;Da_i], or a
%4-row matrix [t_i;a_i;up_Da_i;low_Da_i]; as a 2-column [t_i a_i], a 3-column 
%[t_i a_i Da_i], or a 4-column matrix [t_i a_i up_Da_i low_Da_i]; or as 
%contiguous vectors in the input argument list: relaxfit(...,t_i,a_i,...), 
%relaxfit(...,t_i,a_i,Da_i,...), or relaxfit(...,t_i,a_i,up_Da_i,low_Da_i,...).
%Here t_i, a_i, Da_i are the experimental frequencies and the corresponding 
%amplitudes ad (symmetric) errorbars, respectively; up_Da_i and low_Da_i are 
%the upper and lower (asymmetric) errorbars.
%If the errorbars Da_i are missing, they are loaded from the environment
%variable NMRnoise if its size matches, otherwise they are assumed equal to 1. 
%By default, data are loaded from the environment: 
%[NMRpar2d.val;NMRamp;NMRnoise] (see ENVIRONMENT for details).
%
%
%GUESS    
%Vector, initial guess of fitting parameters. Its length determines the number
%of components of the model. By default, a heuristic initial guess is applied.
%
%
%COMPONENTS
%Scalar, number of additive fitting components (e.g. exponentials). A possible 
%constant offset is not included in the number of components. 
%This parameter is overridden by GUESS. Default is one component.
%
%
%Option SWITCHES:
%Strings with a leading "-" character. Supported options are:
%
%  -w	Force data to be loaded from the workspace environment instead of 
%	the input argument list.
%
%  -b    Force the fit mode of FMINUIT to automatic (bacth). This overrides 
%	the global setup.
%
%  -m    Force the fit mode of FMINUIT to manual (interactive). This overrides 
%	the global setup.
%
%  -a<up_errbar_factor>:<low_errbar_factor>
%  (i.e., "-a" followed by the string representation of two numbers separated by
%   a colon, e.g. '-a1:1.3' )
%	     Use asymmetric error bars. The actual upper and lower error bars 
%	     are obtained as ERRORBARS * up_errbar_factor and ERRORBARS * 
%	     low_errbar_factor, respectively. Here ERRORBARS are the 
%	     symmetric errorbars (either a scalar uniform value or a vector).
%	     If asymmetric errorbars are passed trough either a 4th vector in 
%	     the data vector list, or a 4-row data matrix, or a 2-row NMRnoise 
%	     global matrix, this option switch is ignored.
%	     Asymmetric error bars may be useful when a spectrum shows dips
%	     due to artifacts, in order to prevent the fit to follow such dips.
%
%  -s
%  -c	     Fminuit options. Such switches require an extra option argument
%	     (a matrix or a string, respectively) following in the ZEEMNUQ 
%	     command line. The option switch and the option argument are passed
%	     to fminuit: see fminuit's help for reference. 
%
%
%
%OUTPUT PARAMETERS
%
%BESTPARS      Best-fit parameters
%
%D_BESTPARS    RMS errors of BESTPARS
%
%CHISQUARE     Normalized Chi^2
%
%CORR_MATRIX   Correlation matrix of BESTPARS       
%
%
%
%EXAMPLES
%
%The following is the decriptor for T2FITM (multiexponential decay):
%2 parameters per component (types = time constant, amplitude)+ offset, 
%interlaced parameters.
%   
%>>descr=struct('fcn','t2fitm','np',2,'offs',1,'interl',1,'pnam',...
%              ['T2/2';'Amp '],'pty', 'ta');
%	    	
%
%Fit of a T2 experiment with 2 exponentials (data loaded from the environment):
%
%>> [a b c] = relaxfit(descr,2);
%
%or (an initial guess is specified)
%
%>> [a b c] = relaxfit(descr,[1e-3 5e-3 10 20 0]);
%
%
%Fit of simulated data, inlined through the input list, with 2 exponentials:
%
%>>t=0:100;
%>>y=10*exp(-t/10)+20*exp(-t/30)+.05*(rand(size(t))-.5);
%>>dy=.02*ones(size(t));
%
%>> [a b c] = relaxfit(descr,t,y,dy,2);
%or
%>> [a b c] = relaxfit(descr,[t;y;dy],2);
%
%or (an initial guess is specified)
%
%>> [a b c] = relaxfit(descr,t,y,dy,[4 55 3 4 0]);
%or
%>> [a b c] = relaxfit(descr,[t;y;dy],[4 55 3 4 0]);
%
%
%   NMR2001: relaxfit 
%   Revision: 2.0,  14-Sept-2002
%   Copyright (c) 2001-03 by Giuseppe Allodi
%
%
function [varargout]= relaxfit(rdescr,varargin);
%% RELAXFIT: relaxation fit gateway 

global NMRamp NMRshift NMRnoise;
global NMRpar2d NMRplotsty NMRfitflags
global commands parnames stepbounds;
varargout=cell(1,max(3,nargout)); %=[cf,errs,chi,[emt]]

if (~nargin),
  usedefault=1; 
elseif ~is_rdescr(rdescr), %%if not a descriptor, push arg into vararg list
  varargin=cat(2,{rdescr},varargin);
  usedefault=1;
else
  usedefault=0;
end 
if (usedefault),
  rdescr=struct('fcn','t2fitm','np',2,'offs',1,'interl',1,'pnam',...
    ['T2/2';'Amp '],'pty','ta');
  %%rdescr.np=2; %% No. of parameters for each component
  %%rdescr.offs=1; %% 1: uses offset
  %%rdescr.interl=1; %% 1: interlaced parameters, e.g.: T2#1, T2#2, A#1, A#2
  %%rdescr.pty='ta'; %%'a'=amplitude, 'r'==rate, 't'=time constant, 'b'==beta
end
if ((rdescr.np < 2) | (size(rdescr.pnam,1) < rdescr.np) | ...
    (length(rdescr.pty)<rdescr.np)),
  error('Incomplete parameter description');
end
if isunix,
  fntsz=12;
else
  fntsz=10;
end

%flags
use_wrkspc=0;
batchfit=NMRfitflags.batch;
const_noise=0; %if set, use a constant NMRnoise as uniform errorbars
asymfac=0; %if asymfac=[upfac, lowfac], asymmetric error bars = 
	   %[upfac*bars; lowfac*bars] 
 		 
argc=length(varargin);

if (~argc) %% 1 component, use workspace
  use_wrkspc=1;
end

%%parse inputs
vec_ix=[]; scal_ix=[]; matr_ix=[];
fmn_ix=[]; skipnext=0; 
for k=1:argc;
  if (skipnext),
    skipnext=0;
  elseif isstr(varargin{k}),
    if (strncmp(varargin{k},'-w',2)), %force to get data from workspace
      use_wrkspc=1;
    elseif (strncmp(varargin{k},'-s',2) | strncmp(varargin{k},'-c',2)), 
      %arguments inlined to fminuit
      skipnext=1;
      fmn_ix=[fmn_ix, k:min(k+1,argc)];
    elseif (strncmp(varargin{k},'-b',2)), %batch
      batchfit=1;
      disp('Fit in batch mode')
    elseif (strncmp(varargin{k},'-m',2)), %interactive
      batchfit=0;
      disp('Interactive fit')
    elseif (strncmp(varargin{k},'-a',2)), %error bar asymmetry
      xx=sscanf(varargin{k},'-a%g:%g');
      if (length(xx)==2 & all(xx)),
	asymfac=xx(:);
	disp('Fitting with asymmetric errorbars');
      else
	disp('Invalid ''-a'' option');
      end
    else
      disp(sprintf('Unknown ''%s'' option.',varargin{k}));
    end
  elseif (isscal(varargin{k})), %%the evaluation order is important
    scal_ix=[scal_ix,k];
  elseif (isvec(varargin{k})),
    vec_ix=[vec_ix,k];
  elseif (ismatr(varargin{k})),
    matr_ix=[matr_ix,k];
  end
end  

%%locate data:
if (use_wrkspc | (~length(matr_ix) & (length(vec_ix)<2))), %%% 1) Use workspace
  if ((size(NMRpar2d.val) ~= size(NMRamp)) | ~size(NMRamp,2)),
    error('No data or bad workspace')
  elseif (~(strncmp(NMRpar2d.name2d,'Delay',5) | ...
            strncmp(NMRpar2d.name2d,'Misc',4))),
    error('Not a relaxation experiment');
  end
  if (size(NMRnoise,2) == size(NMRamp,2)), %use NMRnoise as errorbars
    data=[NMRpar2d.val;NMRamp;NMRnoise];
  else
    data=[NMRpar2d.val;NMRamp];
  end
  use_wrkspc=1;
elseif (length(matr_ix) & length(vec_ix)<2),  %%% 2) Data in a matrix
  data=varargin{matr_ix(1)};
  sz=size(data); 
  if (sz(1)>4 & sz(2)<=4), %transpose if data are arranged column-wise
    data=data.';
  elseif (sz(1)>4 & sz(2)>4),
    error('Invalid data matrix');
  end
else  %%% 3) x, y, [Dy], [Dy_low] data in vectors
  len=length(varargin{vec_ix(1)}); 
  if (len ~= length(varargin{vec_ix(2)})),
    error('Size mismatch between X and Y data');
  end
  data=[varargin{vec_ix(1)}(:).';...
	  varargin{vec_ix(2)}(:).'];
  vec_ix(1:2)=[];  %%pop vectors from vararg list

  % error bars 
  if length(matr_ix), %errors in a [Dy_up;Dy_low] matrix are also accepted
    sz=size(varargin{matr_ix(1)});
    if (sz(1)==2 & sz(2)==len),
      data=[data;varargin{matr_ix(1)}];
    elseif (sz(2)==2 & sz(1)==len), %transposed
      data=[data;varargin{matr_ix(1)}'];
    else
      disp('Invalid error matrix');
    end
  end
  for k=1:min(4-size(data,1),length(vec_ix)) %%search errors in vectors
    if (len == length(varargin{vec_ix(1)})),
      data=[data;varargin{vec_ix(1)}(:).'];
      vec_ix(1)=[];
    else
      break;
    end
  end
  
end

sz=size(data); 
if (asymfac), %%asymmetrize errorbars
  if (sz(1)>3) 
    disp('Double errorbars supplied with data, ''-a'' flag ignored');
  elseif (sz(1)==3),
    data=[data(1:2,:); asymfac*data(3,:)];
  else %sz(1)=2
     data=[data; asymfac*ones(1,sz(2))];
     const_noise=1; %flag: no errorbars with data
  end
end


%%look for initial guess or no of component  
n_c=[]; guess=[];
if (length(vec_ix)),
  guess=varargin{vec_ix(1)};
elseif (length(scal_ix)),
  n_c=abs(varargin{scal_ix(1)});
end
if (~(length(n_c) | length(guess))),
  n_c=1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
lef=min(data(1,:)); rig=max(data(1,:));
if (length(guess)),
  n_c=floor(length(guess)/rdescr.np); 
else   %% auto-guess
  %%find occurrence index of amplitudes, rates
  aix=find(rdescr.pty=='a'); laix=length(aix); 
  rix=find(rdescr.pty=='r'); lrix=length(rix);
  tix=find(rdescr.pty=='t'); ltix=length(tix);
  bix=find(rdescr.pty=='b'); lbix=length(bix);
  if (~(laix & (lrix | ltix)) | (laix+lrix+ltix+lbix ~= rdescr.np)),
    error('Invalid parameter type description');
  end
  n_c=floor(n_c);
  guess=zeros(1,rdescr.np*n_c + rdescr.offs);
  %%automatic initial guess
  aa= max(data(2,:))/(laix*n_c);
  [dmy_ ix_]=sort(-[rix tix]); %% pow_of_2 scaled rates
  tt=(rig-lef)./2.0.^ix_(lrix+1:end);
  rr=2.0.^ix_(1:lrix)/(rig-lef);
  for ii=0:n_c-1;
    guess(pindex(ii,n_c,rdescr.np,aix,rdescr.interl)) = aa*ones(1,laix);
    guess(pindex(ii,n_c,rdescr.np,rix,rdescr.interl)) = rr/exp(ii);
    guess(pindex(ii,n_c,rdescr.np,tix,rdescr.interl)) = tt*exp(ii);
    guess(pindex(ii,n_c,rdescr.np,bix,rdescr.interl)) = ones(1,lbix);
  end
end

parnames=''; %% parameter names in fminuit
sufx='';
if (rdescr.interl), %% interlaced parameters
  for jj=1:rdescr.np;
    for ii=1:n_c;
      if (n_c > 1),
	sufx=sprintf('#%d',ii);
      end
      parnames=[parnames,sprintf('%s%s ',strtok(rdescr.pnam(jj,:)),sufx)];
    end
  end
else %%non-interlaced
  for ii=1:n_c;
    if (n_c > 1),
      sufx=sprintf('#%d',ii);
    end
    for jj=1:rdescr.np;
      parnames=[parnames,sprintf('%s%s ',strtok(rdescr.pnam(jj,:)),sufx)];
    end
  end
end  
if (rdescr.offs), 
  parnames=[parnames,' Backgnd']; 
end

if (sz(1)>2 & ~const_noise), %% errorbars with data
  const_noise=1;
elseif (length(NMRnoise)==1), %%use NMRnoise (scalar) as errorbar
  const_noise=NMRnoise;
else   %%no errorbars: use default
  const_noise=1;
  disp(sprintf(['Using errorbar half width=%g. ',...
      'A different value may be assigned to ''NMRnoise'''],const_noise)); 
end

%%fix offset to 0, if flag is set
if (rdescr.offs & ~NMRfitflags.offs),
  stepbound_save_=stepbounds;
  szsb=size(stepbounds);
  if (szsb(2)<2 | or(stepbounds(:,1)>length(guess)) | ...
      or(stepbounds(:,1)<1)),
    stepbounds=[length(guess) 0];
  else
    ix__=find(stepbounds(:,1)==length(guess));
    if (szsb(2)==3),
      stepbounds=[stepbounds(:,1),-ones(szsb(1),1),stepbounds(:,2:3)];
      szsb(2)=4;
    end
    if (length(ix__)),
      stepbounds(ix__(1),2) = 0;
    else
      stepbounds=[stepbounds;[ix__(1) 0 zeros(1,szsb(2)-2)]];
    end
  end
end
      

%%minimize with fminuit
[varargout{:}]=...
  fminuit(rdescr.fcn,'mnplot',guess,data,char(97+batchfit),const_noise,...
	  varargin{fmn_ix}); 
cf=varargout{1};
%%normalize chi^2
norma = size(data,2) - length(cf) + length(find(varargout{2}==0));
varargout{3}=varargout{3}/norma;       %%normalization of chi^2

%%restore stepbounds
if (rdescr.offs & ~NMRfitflags.offs),
  stepbounds=stepbound_save_;
end

%% plot
lef_=max(5e-6,lef); rig_=max(5e-6,rig); 
xx=lef_*exp((0:2.5e-3:1)*log(rig_/lef_));
if (lef < lef_),
  xx=[(lef:(lef_-lef)/250:lef_),xx];
end
yy=feval(rdescr.fcn,cf,xx);


if (NMRplotsty.ebar), %%if  errorbar enabled
  if (sz(1)> 2),
    errorbar(1e6*data(1,:),data(2,:),data(3,:),strtok(NMRplotsty.sym));
  else
    err_=const_noise*ones(1,sz(2));
    errorbar(1e6*data(1,:),data(2,:),err_,strtok(NMRplotsty.sym));
  end
  hold on
  plot(1e6*xx,yy,NMRplotsty.line);
  hold off; 
else
  plot(1e6*xx,yy,NMRplotsty.line, 1e6*data(1,:),data(2,:),NMRplotsty.sym);
end
if (NMRplotsty.ylog),
  set(gca,'yscale','log');
end

%%draw grid, if flag is set
if (NMRplotsty.grid),
  grid;
end

xlabel('Delay [\mus]','fontsize',fntsz);
ylabel('Amplitude [a.u.]','fontsize',fntsz);
title('Relaxation experiment','fontsize',fntsz); 



%%%local macros
function r=is_rdescr(a);
r = isstruct(a);
if (r),
  r = (length(fieldnames(a)) == 6);
end

function r=isscal(a);
r=isnumeric(a) & ~isempty(a) & (max(size(a))==1);

function r=isvec(a);
r=isnumeric(a) & (min(size(a))==1);

function r=ismatr(a);
r=isnumeric(a) & (min(size(a))>1);

function k=pindex(i,nc,np,ix,fl);  %%fl==1: interlaced parameters
if ~length(ix),
  k=[];
  return
end
if (fl),
  k=1+i+(ix-1)*nc;
else
  k=np*i+ix;
end

%% Copyright (C) G. Allodi
%% Revision 2.0   14-sept-2002
