%ZEEMNUQ		Fitting and plot of spectral data to a powder quadrupole 
%		pattern with gaussian broadening for B>>Q. 
%		The quadrupole satellites are calculated up to the 2nd order 
%		in Q/B from perturbation theory (high field approximation). 
%		Powder averaging is performed numerically according to  
%		Alderman's algorithm. Magnetic broadening is simulated by 
%		convolution with a gaussian lineshape.
%		Chi^2 fitting is performed by the FMINUIT minimization engine.
% 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = 
%	      ZEEMNUQ(CONSTANTS,...)
%
%(dots stand for a variable argument list). Optionally, the argument list
%may contain experimental DATA (by default, they are loaded from the 
%environment, see also ENVIRONMENT), initial GUESS of fitting parameters, 
%and option SWITCHES.
%
%
%INPUT PARAMETERS
%
%CONSTANTS
%A 5-element vector, CONSTANTS=[f_min f_width, points, precision, spin],
%where:
%  f_min		lower limit of the frequency window employed for numerical 
%		powder average;
%
%  f_width	width of of the frequency window;
%
%  points	number of points the frequency window is divided into;
%
%  precision	precision parameter; the higher is precision, the more 
%		accurately powder singularities are reproduced, and the 
%		slower is calculation; values between 30 and 60 are usually
%		a good compromise;
%
%  spin		magnitude of the nuclear spin (e.g. 3/2)
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
%An implicit inline method is also implemented by means of the "-z" switch if 
%the Fourier transform of a FID (or spin echo) is to be treated as the 
%spectrum and fitted (see below).
%By default, data are loaded from the environment: 
%[NMRpar2d.val;NMRamp./NMRpar2d.val.^2;NMRnoise./NMRpar2d.val.^2] 
%(see ENVIRONMENT for details). Notice that the amplitude correction by 
%omega^2 is automatically performed in the latter case.
%
%
%GUESS    
%5-element vector, initial guess of fitting parameters. 
%GUESS = [nuL, nuQ, eta, amplitude, sigma], 
%defined as:
%
%nuL	    the Larmor frequency gamma*B/(2*pi);
%
%nuQ	    the quadrupole frequency (the quadrupole coupling parameter: it 
%	    coincides with the zero-field NQR frequency only if eta==0);
%
%eta	    the EFG asimmetry parameter;
%
%amplitude   integral amplitude of the spectrum;
%
%sigma	    magnitude of the gaussian line broadening.
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
%  -r	     Take frequency limits as exact, with no allowed adjustment. This 
%	     is the desired behaviour when fitting a single FID/spin echo 
%	     spectrum. By default, the frequency window limits f_min, f_width
%	     (see above) are internally tuned to match the experimental 
%	     frequency points. The latter is convenient for point-by-point
%	     spectra, where frequency points may not be equally spaced.   
%	     
%  -z[<track>][<zone>] 
%  (i.e., "-z" optionally followed by a character (<track>) and a substring 
%   which is the representation of an integer (<zone>) )
%             The Fourier transform of the FID (or spin echo) at the <zone>th 
%	     zone is inlined as the spectrum. The <track> character may be 
%	     one of "r", "i", or "a", specifying that the real part, imaginary 
%	     part, or modulus of FFT(FID) is to be passed, respectively. 
%	     The frequency base is calculated automatically from the instrument
%	     reference frequency, the FFT analysis shift (stored in the 
%	     NMRshift variable), and the frequency resolution. 
%	     Default is modulus, last zone.
%	     This option implies -r. 
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
%BESTPARS      Best-fit parameters; see above GUESS for parameter definition.
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
%Fit of a spectrum with the  powder quadrupole pattern of a spin 5/2 (data are
%loaded from the environment and corrected for omega^2). The spectrum is 
%calculated over 1024 points, with precision parameter 30:
%
%>> [a b c] = zeemnuq([30 30 1024 30 5/2]);
%
%or (an initial guess is specified)
%
%>> [a b c] = zeemnuq([30 30 1024 30 5/2],[44 3 .2 10 .1]);
%
%
%   NMR2001: zeemnuq 
%   Revision: 1.1,  08-Sept-2002
%   Copyright (c) 2001-02 by Giuseppe Allodi
function [varargout]=zeemnuq(varargin);

global NMRamp NMRshift NMRnoise;
global NMRpar2d NMRplotsty NMRfitflags
global commands parnames stepbounds;

varargout=cell(1,max(3,nargout)); %=[cf,errs,chi,[emt]]

%flags
use_wrkspc=0;
fitFFT=0; 
batchfit=NMRfitflags.batch;
const_noise=0; %if set, use a constant NMRnoise as uniform errorbars
asymfac=0; 
%if asymfac=[upfac, lowfac], asymmetric error bars = [upfac*bars; lowfac*bars] 

if isunix,
  fntsz=12;
else
  fntsz=10;
end

swtchopt='s'; %%default: auto-tune freq interval */
%%parse inputs
ylabstr='Amplitude [a.u.]';
vec_ix=[]; scal_ix=[]; matr_ix=[];
is_veclen_5=[]; %% vector of boolean, true if len(vec) == 5
fmn_ix=[]; skipnext=0; 
for k=1:length(varargin);
  if (skipnext),
    skipnext=0;
  elseif isstr(varargin{k}),
    if (strncmp(varargin{k},'-w',2)),
      use_wrkspc=1;
    elseif (strncmp(varargin{k},'-s',2) | strncmp(varargin{k},'-c',2)), 
      %arguments inlined to fminuit
      skipnext=1;
      fmn_ix=[fmn_ix, k:min(k+1,length(varargin))];
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
      swtchopt='r'; %%raw mode
      global FDfids;
      disp(sprintf('Getting data from %s(FFT) of zone # %d',track,fitFFT));
    elseif (strncmp(varargin{k},'-b',2)), %batch
      batchfit=1;
      disp('Fit in batch mode')
    elseif (strncmp(varargin{k},'-m',2)), %interactive
      batchfit=0;
      disp('Interactive fit')
    elseif (strncmp(varargin{k},'-r',2)),
      swtchopt='r'; %%raw mode
      disp('Disabling frequency window adjustment')
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
    is_veclen_5=[is_veclen_5,(length(varargin{k})==5)];
  elseif (ismatr(varargin{k})),
    matr_ix=[matr_ix,k];
  end
end  


%% locate 5-vector = [freq_win_par,[init. guess]]: keep the rightmost 2 
ix_=find(is_veclen_5); n_=length(ix_); 
ix_=ix_(1:n_-2); %%indices of leftmost 5-vectors: kept as regular vectors
is_veclen_5(ix_)=~ones(1,length(ix_)); %must be a logical op

vec5_ix=vec_ix(is_veclen_5); %% indices of actual special 5-vectors
vec_ix=vec_ix(~is_veclen_5); %%other (regular) vectors
if (~length(vec5_ix)),
  error(['Frequency-window and spin parameters must be supplied in a ',...
      '5-element input vector [f_start,f_width,N-points,N_tents,I].']);
end
setup_pars=varargin{vec5_ix(1)};

%%locate data:
if (length(vec_ix)==1 | length(vec_ix)>3 | (length(matr_ix) & length(vec_ix))),
  warning(' unmatched vector arguments ignored');
end
if (fitFFT), %%%1) Fit fourier transform of a FID
  frq=freqbase(fitFFT);
  if (size(NMRnoise,2) == length(frq)),
    data=[frq;eval([track,'(FDfids(1:length(frq),fitFFT))''']); NMRnoise];
  else
    data=[frq;eval([track,'(FDfids(1:length(frq),fitFFT))'''])];	  
  end
  %In this context, a negative frequency span (in setup_pars) has
  %the special meaning that frequency bounds are calculated from
  %freqbase(fitFFT). The frequency span is doubled.  
  if ((setup_pars(2) < 0) | (setup_pars(3) < 0)),    
    setup_pars(3)=max(2*length(frq),setup_pars(3)); %double N, if needed
    f0=[frq(1)+2*frq(end)-frq(end-1)]/2; %freqbase yields folded shifts
						    
    setup_pars(1)=f0-length(frq)*(frq(1)-frq(2));
    setup_pars(2)=2*length(frq)*(frq(1)-frq(2));
  end
elseif (use_wrkspc | (~length(matr_ix) & (length(vec_ix)<2))), 
  %%% 2) Use workspace
  if (~strcmp(NMRpar2d.name2d,'Frequency')),
    error('Not a frequency-swept spectrum');
  elseif ((size(NMRpar2d.val) ~= size(NMRamp)) | ~size(NMRamp,2)),
    error('No data or bad workspace')
  end
  if (size(NMRnoise,2) == size(NMRamp,2)), %use NMRnoise as errorbars
    data=[NMRpar2d.val; NMRamp; NMRnoise];
  else  %NMRnoise corrupted, use unit conventional errorbars 
    data=[NMRpar2d.val;NMRamp; ones(size(NMRamp))]; 
    const_noise=1; %flag: no errorbars with data
  end
  data(2:end,:)=data(2:end,:) ...
      ./(ones(size(data,1)-1,1)*data(1,:).^2);   %scaling by 1/nu^2 
  use_wrkspc=1;
  ylabstr='Amplitude/\nu^2 [a.u.]';
elseif (length(matr_ix) & length(vec_ix)<2),  %%% 3) Data in a matrix
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
  for k=3:min(6-size(data,1),length(vec_ix)) %%search errors in vectors
    if (len == length(varargin{vec_ix(k)})),
      data=[data;varargin{vec_ix(k)}(:).'];
    else
      warning(' Size mismatch between vector arguments: errorbars ignored');
      break
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

lef=min(data(1,:)); rig=max(data(1,:));
%%look for initial guess 
if (length(vec5_ix)>1),
  guess=varargin{vec5_ix(2)};
else %%auto-guess
  [dummy,ix]=max(data(2,:));
  w=(rig-lef)/6;
  A0 = sum(data(2,:))*mean(diff(sort(data(1,:))));
  guess = [data(1,ix),w,.25,A0,.2];
end

  %%%%%%%%%%%%%%%%%%%%%%%%%%
parnames='nuL nuQ Eta Amp Sigma';

if (sz(1)>2 & ~const_noise), %% errorbars with data
  const_noise=1;
elseif (length(NMRnoise)==1), %%use NMRnoise (scalar) as errorbar
  const_noise=NMRnoise;
else   %%no errorbars: use default
  const_noise=1;
  disp(sprintf(['Using errorbar half width=%g. ',...
      'A different value may be assigned to ''NMRnoise'''],const_noise)); 
end

%%setup nuzqfit function
nuzqfit(swtchopt,setup_pars,data(1,:)); 

%%minimize with fminuit
[varargout{:}]=...
  fminuit('nuzqfit','nuzqplot',guess,data,char(97+batchfit),const_noise,...
	  varargin{fmn_ix}); 
cf=varargout{1};
%%normalize chi^2 
norma = size(data,2) - length(cf) + length(find(varargout{2}==0));
varargout{3}=varargout{3}/norma; %%chi^2 normalization

%%%%%%	plot
[sp,frq]=nuzqfit(cf,1:2); %% best fit curve; 1:2=dummy arg

if (NMRplotsty.ebar),
  if (sz(1)> 2),
    errorbar(data(1,:),data(2,:),data(3,:),data(3,:),NMRplotsty.sym);
  else
    err_=const_noise*ones(1,sz(2));
    errorbar(data(1,:),data(2,:),err_,err_,NMRplotsty.sym);
  end
  hold on
  plot(frq,sp,NMRplotsty.line);
else
  plot(frq,sp,NMRplotsty.line, data(1,:),data(2,:),NMRplotsty.sym);
end
%%draw grid, if flag is set
if (NMRplotsty.grid),
  grid;
end

hold off;
xlabel('Frequency (MHz)','fontsize',fntsz);
%ylabel('Amplitude [a.u.]','fontsize',fntsz);
ylabel(ylabstr,'fontsize',fntsz);
title('NMR Spectrum','fontsize',fntsz);

function r=isscal(a);
r=isnumeric(a) & ~isempty(a) & (max(size(a))==1);

function r=isvec(a);
r=isnumeric(a) & (min(size(a))==1);

function r=ismatr(a);
r=isnumeric(a) & (min(size(a))>1);

%% Copyright (C) G. Allodi
%% Revision 2.0   14-sept-2002
