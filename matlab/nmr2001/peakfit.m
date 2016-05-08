%PEAKFIT		Gateway program for fitting and plot of spectral data. 
%		Spectra are fitted to one or many peaks by the FMINUIT 
%		minimization engine. The multi-component fit model is 
%		specified by a suitable descriptor.  
%		NB.: This routine is not intended to be called directly.
%		A user should call instead the appropriate wrapper function 
%		(e.g. mgau), which takes care of passing the descriptor of 
%		the appropriate fitting function.
% 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = PEAKFIT([DESCRIPTOR],...)
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
%default function is many gaussian peaks.
%
% DESCRIPTOR.fcn (string) 
%  The fit function name. The fit function has the prototype y=fcn(P,D), 
%  where P is a vector of fit parameters, and D is 1-,2-, or 3-row matrix 
%  of experimental data. 
%  The length of P is variable: 
%    length(P) = DESCRIPTOR.np * NC
%  or
%    length(P) = DESCRIPTOR.np * NC + 1
%  where NC is the number of components, and DESCRIPTOR.np is the number of 
%  parameter per component (see below). If length(P)== DESCRIPTOR.np*NC + 1
%  (NC integer), the last element on P is a constant component (offset).
%  The ordering of fit parameters is non-interlaced. In the case for instance
%  of 2-component model, each component depending on 3 free parameters A_,
%  B_, C_, the order is 
%    P = [A_1, B_1, C_1, A_2, B_2, C_2]; 
%  If D is a row-vector, it is interpreted as the frequency-base vector and 
%  fcn() returns the corresponding values of the theoretical function; otherwise
%  it is interpreted as [t;Y_t] or [t;y_t;Delta_Y_t] and the function returns 
%  the unnormalized Chi^2.
%
% DESCRIPTOR.np (scalar) 
%  The number of fit parameters for each component. 
%  For instance, DESCRIPTOR.np=3 (amplitude, center, width) for a multi-
%  gaussian model (GMFIT).  
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
%  'c': peak center;
%  's': peak width;
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
%An implicit inline method is also implemeted by means of the "-z" switch if 
%the Fourier transform of a FID (or spin echo) is to be treated as the 
%spectrum and fitted (see below).
%By default, data are loaded from the environment: 
%[NMRpar2d.val;NMRamp./NMRpar2d.val.^2;NMRnoise./NMRpar2d.val.^2] 
%(see ENVIRONMENT for details). Notice that the amplitude correction by 
%omega^2 is automatically performed in the latter case.
%
%
%GUESS    
%Vector, initial guess of fitting parameters. Its length determines the number
%of components of the model. By default, a heuristic initial guess is applied.
%
%
%COMPONENTS
%Scalar, number of additive fitting components (peaks). A fractional
%component is interpreted as a constant offset (e.g., COMPONENTS = 2.5
%means 2 peaks plus a constant component. This parameter is overridden by GUESS.
%Default is one component.
%
%
%Option SWITCHES:
%Strings with a leading "-" character. Supported options are:
%
%  -w	     Force data to be loaded from the workspace environment instead 
%	     of the input argument list.
%
%  -b	     Force the fit mode of FMINUIT to automatic (bacth). This 
%	     overrides the global setup.
%
%  -m	     Force the fit mode of FMINUIT to manual (interactive). This 
%	     overrides the global setup.
%
%  -z[<track>][<zone>] 
%  (i.e., "-z" optionally followed by a character (<track>) and a substring 
%  which is the representation of an integer (<zone>) )
%             The Fourier transform of the FID (or spin echo) at the <zone>th 
%	     zone is inlined as the spectrum. The <track> character may be 
%	     one of "r", "i", or "a", specifying that the real part, imaginary 
%	     part, or modulus of FFT(FID) is to be passed, respectively. 
%	     The frequency base is calculated automatically from the reference 
%	     frequency, the analysis shift (in the NMRshift variable),
%	     and the frequency resolution. 
%	     Default is modulus, last zone.
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
%The following is the decriptor for GMFIT (many gaussian peaks):
%3 parameters per component, of type  amplitude, center, width
%   
%>>descr=struct('fcn','gmfit','np',3,'pnam',['Amp  ';'Cent ';'Sigma'],'pty', ...
%              'acs');
%
%Fit of a spectrum with 2 gaussian peaks (data loaded from the environment 
%and corrected for omega^2):
%
%>> [a b c] = peakfit(descr,2);
%
%or (an initial guess is specified)
%
%>> [a b c] = paekfit(descr,[10 30 2 20 50 1.5]);
%
%Fit of a spectrum with 2 gaussian peaks plus a constant offset (data loaded 
%from the environment and corrected for omega^2):
%
%>> [a b c] = peakfit(descr,2.5);
%
%or (an initial guess is specified)
%
%>> [a b c] = paekfit(descr,[10 30 2 20 50 1.5 .03]);
%
%
%Fit of simulated data, inlined through the input list, with 2 gaussians:
%
%>>x=0:100;
%>>y=10*exp(-.5*((x-30)/2).^2) + 20*exp(-.5*((x-30)/2).^2) + ...
%>>			      .05*(rand(size(x))-.5);
%>>dy=.02*ones(size(t));
%
%>> [a b c] = peakfit(descr,x,y,dy,2);
%or
%>> [a b c] = peakfit(descr,[x;y;dy],2);
%
%or (an initial guess is specified)
%
%>> [a b c] = peakfit(descr,x,y,dy,[4 55 3 4 35 2]);
%or
%>> [a b c] = peakfit(descr,[x;y;dy],[4 55 3 4 35 2]);
%
%
%   NMR2001: peakfit 
%   Revision: 2.2,  22-Feb-2006
%   Copyright (c) 2001-06 by Giuseppe Allodi
%
%
function [varargout]= peakfit(pdescr,varargin);
%% multi-peak fit (default: 1 gauassian)
%% peak definition is entered through the structure 'pdescr'
global NMRamp NMRnoise FDfids;
global NMRpar2d NMRplotsty NMRfitflags
global commands parnames stepbounds;
varargout=cell(1,max(3,nargout)); %=[cf,errs,chi,[emt]]

if (~nargin),
  usedefault=1;
elseif ~is_pdescr(pdescr), %%if not a descriptor, push arg into vararg list
  varargin=cat(2,pdescr,varargin);
  usedefault=1;
else
  usedefault=0;
end  
if (usedefault),
  pdescr=struct('fcn','gmfit','np',3,'pnam',['Amp  ';'Cent ';'Sigma'],...
		'pty','acs');
end
if ((pdescr.np < 3) | (size(pdescr.pnam,1) < pdescr.np) | ...
    (length(pdescr.pty)<pdescr.np)),
  error('Incomplete parameter description');
end

%flags
use_wrkspc=0;
fitFFT=0; 
batchfit=NMRfitflags.batch;
const_noise=0; %if set, use a constant NMRnoise as uniform errorbars
asymfac=0; %if asymfac=[upfac, lowfac], asymmetric error bars = 
           %[upfac*bars; lowfac*bars] 
		 
argc=length(varargin);
if isunix,
  fntsz=12;
else
  fntsz=10;
end

if (~argc) %% 1 component, use workspace
  use_wrkspc=1;
end
%%parse inputs
xlabstr='Frequency [MHz]';
ylabstr='Amplitude [arb.u.]';
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
    elseif (strncmp(varargin{k},'-z',2)), %Fit a FFT track
      fitFFT=sscanf([int2str(length(NMRamp)),' ',...
		  strtok(varargin{k},'-zria')],'%g ');
      fitFFT=fitFFT(end);
      track='abs';
      if (any(find(varargin{k}=='r'))),
	track='real';
      elseif (any(find(varargin{k}=='i'))),
	track='imag';
      end
      disp(sprintf('Getting data from %s(FFT) of zone # %d',track,fitFFT));
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
      disp(sprintf('Unknown ''%s'' option',varargin{k}));
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
if (fitFFT), %%% 1) Fit fourier transform of a FID
  frq=freqbase(fitFFT);
  if (size(NMRnoise,2) == length(frq)),
    data=[frq;eval([track,'(FDfids(1:length(frq),fitFFT))''']); NMRnoise];
  else
    data=[frq;eval([track,'(FDfids(1:length(frq),fitFFT))'''])];	  
  end
elseif (use_wrkspc | ...
	(~length(matr_ix) & (length(vec_ix)<2))), %%% 2) use workspace
  if ((size(NMRpar2d.val) ~= size(NMRamp)) | ~size(NMRamp,2)),
    error('No data or bad workspace')
  elseif (strncmp(NMRpar2d.name2d,'Delay',5)),
    error('Not a frequency-swept spectrum');
  end
  if (size(NMRnoise,2) == size(NMRamp,2)), %use NMRnoise as errorbars
    data=[NMRpar2d.val; NMRamp; NMRnoise];
  else  %NMRnoise corrupted, use unit conventional errorbars 
    data=[NMRpar2d.val;NMRamp; ones(size(NMRamp))]; 
    const_noise=1; %flag: no errorbars with data
  end

  if (strcmp(NMRpar2d.name2d,'Frequency')),
    data(2:end,:)=data(2:end,:) ...
	./(ones(size(data,1)-1,1)*data(1,:).^2);   %scaling by 1/nu^2 
    ylabstr='Amplitude/\nu^2 [arb.u.]';
  else
    disp(sprintf(['Abscissa type "%s" doesn''t warrant it''s frequency: '...
		  'using unscaled amplitude.'], NMRpar2d.name2d));
    xlabstr=['Miscellaneous parameter', NMRpar2d.name2d(5:end)];
  end
  use_wrkspc=1;
elseif (length(matr_ix) & length(vec_ix)<2),  %%% 3) data in a matrix
  data=varargin{matr_ix(1)};
  sz=size(data); 
  if (sz(1)>4 & sz(2)<=4), %transpose if data are arranged column-wise
    data=data.';
  elseif (sz(1)>4 & sz(2)>4),
    error('Invalid data matrix');
  end
else  %%% 4) x, y, [Dy], [Dy_low] data in vectors 
  len=length(varargin{vec_ix(1)}); 
  if (len ~= length(varargin{vec_ix(2)})),
    error('Size mismatch between X and Y data');
  end
  data=[varargin{vec_ix(1)}(:).';...
	  varargin{vec_ix(2)}(:).'];
  vec_ix(1:2)=[];  %%pop vectors from vararg list
  %error bars 
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
  n_c=floor(length(guess)/pdescr.np); 
  use_offs=(rem(length(guess),pdescr.np)>0);
else
  aix=findstr(pdescr.pty,'a'); laix=length(aix); %amp. params
  cix=findstr(pdescr.pty,'c'); lcix=length(cix); %position params
  six=findstr(pdescr.pty,'s'); lsix=length(six); %width params
  bix=findstr(pdescr.pty,'b'); lbix=length(bix); %adimensional params
  if (~(laix & lcix & lsix) | (laix+lcix+lsix+bix < pdescr.np)),
    error('Invalid parameter type description');
  end
  use_offs=(rem(n_c,1)>0); n_c=floor(n_c);
  guess=zeros(1,pdescr.np*n_c + use_offs);
  %%automatic initial guess
  ss=(rig-lef)/(4*n_c);
  aa= max(data(2,:))*sqrt(2*pi/n_c)*ss/laix;
  for ii=0:n_c-1;
    cc=lef+ (ii+1)*(rig-lef)/(n_c +1);
    guess(pdescr.np*ii+aix) = aa*ones(1,laix);
    guess(pdescr.np*ii+cix) = cc*ones(1,lcix);
    guess(pdescr.np*ii+six) = ss./2.0.^(0:lsix-1);
    guess(pdescr.np*ii+bix) = ones(1,lbix);
  end
end

parnames=''; %% parameter names in fminuit
for ii=1:n_c;
  for jj=1:pdescr.np;
    parnames=[parnames,sprintf('%s#%d ',strtok(pdescr.pnam(jj,:)),ii)];
  end
end  
if (use_offs), 
  parnames=[parnames,' Backgnd']; 
end

if (sz(1)>2 & ~const_noise), %% errorbars with data
  const_noise=1;
elseif (length(NMRnoise)==1), %%use NMRnoise (scalar) as errorbar
  const_noise=NMRnoise;
else   %%no errorbars: use default
  const_noise=1;
  disp(sprintf(['Using errorbar half width=%g. A different value may be' ...
		' assigned to ''NMRnoise'''],const_noise));
end
  
%%minimize with fminuit
[varargout{:}]=...
  fminuit(pdescr.fcn,'mnplot',guess,data,char(97+batchfit),const_noise,...
	  varargin{fmn_ix}); 
cf=varargout{1};
%%normalize chi^2
norma = size(data,2) - length(cf) + length(find(varargout{2}==0));
varargout{3}=varargout{3}/norma; %%chi^2 normalization

%% plot
n=max(500,size(data,2));
xx=lef:(rig-lef)/n:rig; 
yy=feval(pdescr.fcn,cf,xx);

if (NMRplotsty.ebar),
  if (sz(1)> 2),
    errorbar(data(1,:),data(2,:),data(3,:),strtok(NMRplotsty.sym));
  else
    err_=const_noise*ones(1,sz(2));
    errorbar(data(1,:),data(2,:),err_,strtok(NMRplotsty.sym));
  end
  hold on
  plot(xx,yy,NMRplotsty.line);
else
  plot(xx,yy,NMRplotsty.line, data(1,:),data(2,:),NMRplotsty.sym);
end
hold on;

%%draw fitting component, if flag is set
if (NMRplotsty.drawpeaks & ((n_c > 1) | use_offs)), 
  for ii =0:(n_c-1);	
    y_=feval(pdescr.fcn,cf(pdescr.np*ii+(1:pdescr.np)),xx);
    plot(xx,y_,'-.');
  end
  if (use_offs), 
    plot(xx([1,end]),ones(1,2)*cf(pdescr.np*n_c+1),'--');
  end
end
if (NMRplotsty.grid),
  grid;
end

hold off;
xlabel(xlabstr,'fontsize',fntsz);
ylabel(ylabstr,'fontsize',fntsz);
title('NMR Spectrum','fontsize',fntsz); 


%%%local macros
function r=is_pdescr(a);
r = isstruct(a);
if (r),
  r = (length(fieldnames(a)) == 4);
end

function r=isscal(a);
r=isnumeric(a) & ~isempty(a) & (max(size(a))==1);

function r=isvec(a);
r=isnumeric(a) & (min(size(a))==1);

function r=ismatr(a);
r=isnumeric(a) & (min(size(a))>1);

%% Copyright (C) 2001-06 G. Allodi
%% Revision 2.2   22-feb-2006
