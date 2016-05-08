%ZNUQMBIGAU	Fitting and plot of spectral data to a powder quadrupole 
%		pattern with gaussian broadening for B>>Q, plus many bigaussian
%		components. By bigaussian it is meant the superposition
%		of two concentric gaussian curves (each one depends therefore 
%		on 5 parameters). The quadrupole pattern is calculated in the 
%		same way as in ZEEMNUQ. 
%		Chi^2 fitting is performed by the FMINUIT minimization engine.
%		This specialized function is intended for highly inhomogeneous 
%		quadrupole spectra; the bigaussian shape mimics a quadrupole
%		pattern with brodening sigma comparable to nuQ.
% 
%Usage:
%
%>> [BESTPARS, D_BESTPARS, CHISQUARE, CORR_MATRIX] = 
%	      ZNUQMBIGAU(CONSTANTS,...)
%
%(dots stand for a variable argument list). Optionally, the argument list
%may contain experimental DATA, initial GUESS of fitting parameters or 
%number of COMPONENTS, and option SWITCHES.
%
% 
%INPUT PARAMETERS
%The same as of ZEEMNUQ, except for the following: 
%
%GUESS    
% (5+N*5)-element vector, where N is a non-negative integer: initial guess of 
%  fitting parameters. 
%  GUESS = [nuL, nuQ, eta, amplitude, sigma, Ag1#1, f0#1, ..., 
%	 Ag1#N, f0#N, Sg1#N, Ag2#N, f0#N, Sg2#N], 
%  with:
%  nuL, nuQ, eta, amplitude, sigma: parameters of the quadupole powder pattern,
%   defined as in ZEEMNUQ;
%  Ag1#k:   Integral amplitude of the 1st gaussian of the k-th bigaussian 
%	   component; 
%  f0#k:    Center of the k-th bigaussian component; 
%  Sg1#k:   Width of the 1st gaussian of the k-th bigaussian component; 
%  Ag2#k:   Integral amplitude of the 2nd gaussian of the k-th bigaussian 
%	   component; 
%  Sg2#k:   Width of the 2nd gaussian of the k-th bigaussian component; 
%
%COMPONENTS
%Scalar, number of additive gaussian components in addition to the quadrupole 
%powder pattern. This parameter is overridden by GUESS.
%Default is one component.
%
%OUTPUT PARAMETERS
%The same as of ZEEMNUQ.
%
%
%   NMR2001: znuqmbigau 
%   Revision: 1.1,  03-Feb-2003
%   Copyright (c) 2001-02 by Giuseppe Allodi
function [varargout]=znuqmbigau(varargin);

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
asym_ebar=1; %default:symmetric error bars
is_veclen_5=[]; %% vector of boolean, true if len(vec) == 5
is_veclen_10=[]; %% vector of boolean, true if len(vec) == 5+5*n, n>0
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
    is_veclen_10=[is_veclen_10,length(varargin{k})>5 & ...
			      ~rem(length(varargin{k}),5)];
  elseif (ismatr(varargin{k})),
    matr_ix=[matr_ix,k];
  end
end  


%% locate 5-vector = [freq_win_par,[init. guess]]: keep the rightmost 2 
ix_=find(is_veclen_5); n_=length(ix_); 
ix_=ix_(1:n_-2); %%indices of leftmost 5-vectors: kept as regular vectors
is_veclen_5(ix_)=~ones(1,length(ix_));  %must be a logical op

vec5_ix=vec_ix(is_veclen_5); %% indices of actual special 5-vectors

if (~length(vec5_ix)),
  error(['Frequency-window and spin parameters must be supplied in a ',...
      '5-element input vector [f_start,f_width,N-points,N_tents,I].']);
end
setup_pars=varargin{vec5_ix(1)};

%%look for initial guess 
guess=[];
vec10_ix=vec_ix(is_veclen_10); %% indices of (5+5*n)-vectors

if (length(vec5_ix)>1),
  guess=varargin{vec5_ix(2)};
  is_veclen_10=(is_veclen_10 & 0);  %discard all 10-vectors
else %look for guess among (5+5*n)-vectors
  if ((length(vec10_ix) == 1) | (length(vec10_ix) > 3)),
    guess=varargin{vec10_ix(1)};  
  elseif (length(vec10_ix) > 1),
    if (length(varargin{vec10_ix(1)}) ~= length(varargin{vec10_ix(2)})),
      guess=varargin{vec10_ix(1)};  
    end
  end
  
  %discard all 10-vectors but the initial guess   
  ix_=find(is_veclen_10); n_=length(ix_); 
  ix_=ix_(1+(length(guess)>0):n_); %%indices of discarded 10-vectors
  is_veclen_10(ix_)=~ones(1,length(ix_)); %must be a logical op
end

%%keep only other (regular) vectors
vec_ix=vec_ix(~(is_veclen_5 | is_veclen_10)); 

%%locate data:
if (length(vec_ix)==1 | length(vec_ix)>3 | (length(matr_ix) & length(vec_ix))),
  warning(' unmatched vector arguments ignored');
end
if (fitFFT), %%% 1) Fit fourier transform of a FID
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
elseif (use_wrkspc | ...
	(~length(matr_ix) & (length(vec_ix)<2))), %%% 2) Use workspace
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
elseif (length(matr_ix) & length(vec_ix)<2),  %%% 3) Data in a matrix),  
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
%%look again for initial guess 
if (~length(guess)), %no guess inlined, enable autoguess 
  nbgc = 1; %default: 1 gaussian 
  if (length(scal_ix))
    nbgc=abs(round(varargin{scal_ix(1)}));
  end
  [dummy,ix]=max(data(2,:));
  w=(rig-lef);
  A0 = sum(data(2,:))*mean(diff(sort(data(1,:))))/(nbgc+1);
  guess = [data(1,ix),w/6,.25,A0,.1];
  for k=1:nbgc;
    guess=[guess,A0/2,lef+k*w/(nbgc+1),.25*w/(nbgc+1),A0/2,.5*w/(nbgc+1)];
  end
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

for ii =1:(length(guess)-5)/5;
  parnames=[parnames,sprintf(' A1g#%d f0#%d S1g#%d A2g#%d S2g#%d',ii,ii,...
			     ii,ii,ii)];
end
  %%setup nuzqmbgfit function
nuzqmbgfit(swtchopt,setup_pars,data(1,:)); 

%%minimize with fminuit
[varargout{:}]=...
  fminuit('nuzqmbgfit','nuzqplot',guess,data,char(97+batchfit),const_noise,...
	  varargin{fmn_ix}); 
cf=varargout{1};
%%normalize chi^2 
norma = size(data,2) - length(cf) + length(find(varargout{2}==0));
varargout{3}=varargout{3}/norma; %%chi^2 normalization

%%%%%%	plot
[sp,frq]=nuzqmbgfit(cf,1:2); %% best fit curve; 1:2=dummy arg
sz=size(data); 

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

%%draw fitting component, if flag is set
if (NMRplotsty.drawpeaks & (length(cf) > 5)),
  hold on;
  [sp,frq]=nuzqmbgfit(cf(1:5),1:2); %% quadrupolar component; 1:2=dummy arg
  plot(frq,sp,'-.');
  
  for ii =0:(length(cf)-5)/5-1;	
    y_=mbgfit(cf(5+(1:5)+ii*5),frq);
    plot(frq,y_,'--');
  end
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
%End main function

%Local macros
function r=isscal(a);
r=isnumeric(a) & ~isempty(a) & (max(size(a))==1);

function r=isvec(a);
r=isnumeric(a) & (min(size(a))==1);

function r=ismatr(a);
r=isnumeric(a) & (min(size(a))>1);

%% Copyright (C) G. Allodi
%% Revision 2.0   14-sept-2002
