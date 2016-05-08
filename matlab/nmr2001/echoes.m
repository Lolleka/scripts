%ECHOES		Analysis of NMR data: apodization and Fourier transform. 
%		Supported formats are TecMag, Stelar (all versions) and Parma. 
%
%Usage: 
% 
%>> ECHOES(FILENAME,...)
%
%(dots stand for a variable argument list). The argument list may contain 
%several OPTIONS switches.
%
%
%INPUT PARAMETERS (all input arguments are strings)
%
%FILENAME the name of the raw data file(s), passed to the appropriate data 
% loader (tmgload or  slrload). 
% For TecMag binary files, multiple file search is supported, i.e. FILENAME
% may contain wildcards. See TMGLOAD for details.
% 
%OPTIONS have the general syntax '-<SWITCH>[<SUBOPTIONS>]', where SWITCH 
% is a character (e.g. -p). The SUBOPTIONS substring have the format 
% '[c]<param1>[:<param2>]' (no space inserted), where <param1>, <param2> are 
% separated by ':' and represent numeric parameters (usually an interval). 
% If <param2> is omitted, suitable defaults are applied. 
% A leading 'c' character (when applicable) signifies that parameters are 
% measured in an alternative unit (usually, channel units instead of seconds 
% or Hz). 
% If suboptions are not specified, the corresponding parameters are selected
% interactively via the Windows/XWindow interface.
%
% OPTION OVERVIEW 
% Optional parameters are included by square brackets 
%
%
% -M[<instrument_flag>]   
%   Data type: -M0 -> TecMag .tnt data (default), -M1 -> Stelar (all versions),
%              -M2 -> Parma .pna archive 
%
% -f
%   Produce frequency-domain spectra (apodization and Fourier transform) and 
%   amplitude evaluation from Fourier-transformed data. 
%   If '-f' is not accompanied by the '-t' switch (see below), previous 
%   time-domain data are re-analyzed, i.e. the data loader subroutine is 
%   skipped (the <FILENAME> argument is ignored in that case). Default is
%   -t -f.
%   
%
% -F[[c]<shift>[:<width>]]  (default []: interactive)
%   Determines how to assign an amplitude to each Fourier-transformed FID 
%   (stored in the global matrix FDfids) in a 2-D experiment. Amplitudes are 
%   stored in the global vector MNRamp. The amplitude is taken as the maximum 
%   of abs(FDfids) ( real(FDfids) with the -h switch, see below) in the 
%   frequency interval    [<shift>-<width> <shift>+<width>].
%   If a <shift> is specified but <width> is omitted, the latter is got from 
%   the global struct NMRfftwin.
%   If no sub-option is supplied, the frequency interval is selected 
%   interactively.
%   The frequency shift of the selected point in the Fourier-transformed FID 
%   is returned in the global vector NMRshift.
%   If <shift> is non-zero, the Fourier-transformed FIDs are shifted (global 
%   matrix FDfids is affected). A fractional shift is interpolated. The 'c' 
%   suboption specifies channel units instead of Hz.
%   NB: without the -F switch the frequency interval is defined by the 
%   global struct NMRfftwin.
%   EXAMPLES:
%   -F	     Interactive selection, the maximum of the shaded part of the 
%	     curve is chosen.
%
%   -F1e4:3e4  Maximum is chosen in the [-20 40] kHz window; spectra in 
%	      TDfids are shifted by 10 kHz.
%
%   -Fc1:6     Maximum is chosen between channels -5 and + 7; spectra in 
%	      TDfids are shifted by 1 channel.
%
%
% -h[[c]<phase1>[:<phase2>]]  (default []: interactive)
%   Rephasing of spectra after Fourier transform: Fourier transformed FIDs 
%   are multiplied by 
%   exp(i*(<phase1> + (channel-central_channel) * <phase2>)).
%   The amplitude of each 2-D point is taken form real(FDfids) after 
%   rephasing, instead of abs(FDfids).   
%   If <phase2> is omitted, it is set to 0 by default.
%   The 'c' suboption specifies that phasing coefficients are in degrees 
%   (degrees/channel) instead of radiants (radiants/channel).
%   If no sub-option is supplied, spectra are rephased interactively.
%   EXAMPLES:
%   -h	      Manual (interactive) phasing, amplitudes taken from 
%	      real(FDfids).
%	      NB: Manual phasing also allows distinct phasing 
%	      coefficients for each 2-D point. 
%
%   -h0.7:.002 Phasing coefficients = 0.7 rad, 2e-3 rad/channel for all 2-D 
%	      points; amplitudes taken from real(FDfids). 
%
%   -hc15      Phasing coefficients = 15 degrees, 0 degrees/channel for all 
%	      2-D points; amplitudes taken from real(FDfids). 
%
%
% -t
%   Force loading of time-domain data. If '-t' is not accompanied by the 
%   '-f' switch (see above), subsequent processing is skipped: matrices 
%   containing frequency-domain data are blanked, and the program exits.
%   Default behavior is -t -f.
%   
%
% -T[[c]<position>[:<width>]]   (default []: interactive)
%   Places an apodization window centered at <position> from beginning of 
%   the time-domain FID, whose characteristic width in <width>. Default units 
%   are seconds; if the 'c' suboption is supplied, units are channels.
%   Apodization window properties (window type, default width, etc.) are 
%   controlled by the global struct NMRapod. 
%   If no suboption is supplied, apodization window is selected interactively 
%   (this also enables baseline correction and pre-blanking, see -b and -k 
%   below).
%   NB: without the -T switch, a default apodization window is positioned 
%   automatically, based on the pulse sequence. That algorithm however may be
%   inaccurate, and using -T in all cases is recommended. 	       
%   EXAMPLES:
%   -T	       Interactive apodization. The window is placed by cursor #2,
%	       pre-blanking range by cursor #1, and baseline evaluation 
%	       range by cursor #3.
%
%   -T2e-5:3e-5 Centers a 3e-5-second wide window at 2e-5 seconds from 
%	       beginning of FID.	     
%
%   -T2e-5      Centers a default apodization window at 2e-5 seconds from 
%	       beginning of FID. Defaults are defined by the global struct
%	       NMRapod.  	     
%
%   -Tc20:30    Places a 30-channel wide apodization window centered on 
%	       channel 20 of FID.
%
% -b[c]<left_limit>[:<right_limit>]   
%   Baseline correction. Baseline is evaluated between <left_limit> and 
%   <right_limit>, and subtracted from time-domain fids before Fourier 
%   transform. This operation however is performed on a temporary buffer, 
%   i.e. time-domain data stored in TDfids are not affected.
%   If <right_limit> is omitted, its default value is +Inf (i.e., baseline 
%   is estimated from all channels following <left_limit>)  
%   A 'c' character switches to channel units instead of seconds. 
%   EXAMPLES:
%   -b5e-5:7e-5 Baseline is estimated from 5e-5 to 7e-5 seconds following
%	       the beginning of FID acquisition.  
%
%   -b5e-5      Baseline is estimated from 5e-5 seconds following the 
%	       beginning of FID acquisition, to the end of FID track.  
%
%   -bc100      Baseline is estimated from channel 100 to the end of FID track.
%
%
% -k[c]<limit>
%   FID pre-blanking. This option is aimed at the suppression of spurious 
%   spikes (ringing) from the initial stretch of FIDs.
%   Channels at earlier times than <limit> are zeroed before Fourier transform. 
%   This operation is performed on a temporary buffer and does not affect 
%   time-domain FIDs stored in TDfids.
%   if the 'c' suboption is present, <limit> is expessed in channels. 
%   EXAMPLES:
%   -k5e-6      Earlier channels than 5e-6 seconds from the beginning of 
%	       FID acquisition are blanked. 
%
%  -kc10        Channels 1:9 are blanked.
%
%
% -z[<FFTsize>]
%   Zero padding. The buffer for FFT calculation is inflated to 1:FFTsize
%   (time domain FIDs in TDfids are unaffected), with trailing channels set
%   to 0. If FID size is <td_size>, spectral resolution is artificially 
%   increases by <FFTsize>/<td_size>.
%   EXAMPLES:
%   -z1024      If FID length is N < 1024, 1024-N zeroes are appended to 
%	       FID before Fourier transform. 
%
%   -z	       This is equivalent to -z<Max_FID_length>, where 
%	       <Max_FID_length> is the maximum length of FIDs (FID sizes may
%	       be different). If all FIDs have the same length, this has 
%	       no effect. 
%     	       
%
% -I
%   Incremental loading. New data are loaded, analyzed and appended to a 
%   previous dataset, instead of replacing it. 
%
%
% -S
%   Stimulated echo. This option only affects default positioning of the 
%   apodization window. The latter is overridden by the -T option.
%
% -n
%   Noise level estimate (in the Fourier transform domain). Noise is estimated 
%   from the noisy tail of FIDs. The usefulness of this option is dubious; it 
%   might yield a sensible estimate for the uncertainty of 2-D point amplitudes
%   (NMRamp) only for 2-D experiments at fixed frequency.
%   If -b is supplied, noise is estimated over the same interval employed for
%   baseline evaluation; otherwise it is estimated over a default interval,
%   namely [2/3 1] of FID time window. 
%   The noise level is stored in the global vector NMRnoise.
%   Without -n, a conventional value divided by sqrt(N_scans) is returned 
%   to NMRnoise.
%
%
% -p
%   Apply a previous analysis (by employing the same options for apodizations, 
%   etc.) to the current dataset. 
%   NB: the order of -p in the option sequence is relevant. -p overrides 
%   preceeding option switches, whereas it is overridden by option switches
%   that follow it.
% 
%
% machine-specific OPTIONS: 
% -A<acronym>  (Stelar, Parma, Ac300) 
% -s<zone_No>  (Stelar, Parma) 
% -e<zone_No>  (Stelar, Parma) 
% -u           (Stelar, Parma) 
% -m           (Stelar)
% -r           (Parma) 
% -q           (Parma)
% In Stelar and Parma mode, these options are passed to the data loader 
%subroutine. See SLRLOAD and PNALOAD for details.
%
%
%OUTPUT ARGUMENTS 
% None explicit. Result of analysis is returned via global matrices. 
% See ENVIRONMENT for a detailed description.
%
%   NMR2001: echoes 
%   Revision: 2.2,  22-Feb-2006
%   Copyright (c) 2001-06 by Giuseppe Allodi
%
function echoes(fname,varargin);
%%echo analysis of Tecmag (.tnt), Stelar (all formats), Parma (.pna) datasets
%% -M0, -M1, -M2  selects Tecamag, Stelar, or Parma,  respectively.
%% Last revised: 20 MAY 04
global FDfids TDfids NMRseq NMRpath;
global NMRamp NMRapod NMRfftwin NMRstat;
global NMRshift NMRnoise NMRpar2d NMRmisc;

%option switches: -

if (~nargin),
  fname='-';
end
if (strncmp(fname,'-',1)),
  disp('Usage:');
  disp('   echoes <filename_pattern> [<options>]'); 
  disp('')
  error(' ');
end
Intrument=0; %%instrument: default=0 (tecmag)
slropts='-A -e -s -u -m';
pnaopts='-A -e -s -u -r -q';
mach_optlist=cell(0); %%-A -e -s -u -m are dedicated Stelar options

%%option flags: initilized 
t_fl=0; f_fl=0; b_fl=0;  %%TD only; FFT from old TD data; baseline correction; 
k_fl=0; S_fl=0; n_fl=0;  %%echo preblanking; SSE; noise estimate;
z_fl=0; Phcf=0; %%zero padding; no FFT phasing 
T_fl=[0 0]; %%T_fl(1): interactive apod. win; T_fl(2): units (1=time, 2=chann.)
F_fl=[0 (1 + NMRfftwin.funit)]; %%F_fl(1): inter. FD win sel.; F_fl(2): units
I_fl=0; %%p_fl=0;  %%incremental data loading; repeat previous analysis 
b_lim=[]; pblank=[];
%%defaults
Twin=[];
Fwin=[NMRfftwin.foffs,NMRfftwin.fwid]; 
apod_p=NMRapod; %%local apodization params. 

%% Process arguments
for ii=1:length(varargin)
  ar=varargin{ii}; l_ar=length(ar);
  if (~isstr(ar)), 
    error('All input arguments must be strings');
  end
  opt = [ar,'  ']; opt=opt(1:2);
  %% Check options
  if (opt == '-M'),  %%instrument: -M0=Tecmag, -M1=stelar, -M2=Parma
    Instrument = round(getopt(ar,1));
    if (Instrument < 0 | Instrument > 3),
      error(sprintf('Invalid intrument type #%d', Instrument)); 
    end
  elseif (opt == '-f'), %%FFT from old TD data
    f_fl=1;
  elseif (opt == '-F'), %% override FFT window
    [Fwin fl_] = getopt(ar,2);
    if (length(Fwin)),
      F_fl(2)=1+fl_; 
      if (length(Fwin) == 1),
        Fwin = [Fwin 0];
      end
    else
      F_fl(1) = 1; %% choose FFT win interactively through Xwin interface 
    end  
  elseif (opt == '-h'), %% rephasing, get NMRamp from real(FDfids)
    [Phcf fl_] = getopt(ar,2); %-hc<x>:<y>: degs; -h<x>:<y> : rads
    if (length(Phcf)),
      if (length(Phcf) < 2),
        Phcf = [Phcf 0];
      end;
      if (fl_),
        Phcf = Phcf*pi/180; %deg 2 rad
      end
    else
      Phcf = -1; %length(Phcf)==1 & Phcf~=0: interactive phasing
    end
  elseif (opt == '-t'), %%TD only
    t_fl=1;
  elseif (opt == '-T'), %% place apodization window
    [Twin fl_] = getopt(ar,2);
    if (length(Twin)),
      T_fl(2)=1+fl_;
    else
      T_fl(1)=1; %% choose apodizationt window interactively 
    end
  elseif (opt == '-b'), %%baseline correction
    [b_lim fl_] = getopt(ar,2);
    if (length(b_lim)),
      b_fl=1+fl_;
      if (length(b_fl) < 2),
        b_lim=[b_lim inf];
      end
    else
      disp(sprintf('Invalid %s option, ignored.',ar)); 
    end 
  elseif (opt == '-k'), %%echo preblanking
    [pblank fl_] = getopt(ar,1);
    if (length(pblank)),
      k_fl=1+fl_;
    else
      disp(sprintf('Invalid %s option, ignored.',ar)); 
    end
  elseif (opt == '-z'), %%zero padding: if no value given, use the full TDsids
    v_ = getopt(ar,1);
    if (length(v_)),
      z_fl=v_;
    else
      z_fl=-1;
    end
  elseif (opt == '-I'), %%increamental loading
    I_fl=~isempty(NMRseq); %False if no dataset is loaded   
  elseif (opt == '-S'), %%SSE (this only affects default position of apod window)
    S_fl=1; 
  elseif (opt == '-n'), %% estimate noise from the noisy fid tail
    n_fl=1;          
  elseif (opt == '-p'), %% keep previous analysis method
    %%NOTE: for -p, the order is important. Options following -p override previous setting         
    if (~iscell(NMRstat) | (length(NMRstat)<12)),
      error('No previous analysis');
    else
      [b_fl,b_lim,k_fl,pblank,n_fl,z_fl, T_fl,Twin, F_fl, Fwin, apod_p, Phcf] = ...
        deal(NMRstat{:});
    end
    disp('Applying the previous analysis to the present dataset') 
  elseif (length(findstr(slropts,opt)) ...
	  | length(findstr(pnaopts,opt))), %%specialized Stelar/Parma option
    mach_optlist{end+1} = ar;
  else
    disp(sprintf('Unknown ''%s'' option. Ignored',opt));
  end
end %%End FOR loop

if (~Instrument), %%TecMag
  if (length(mach_optlist)),
    msg='meaningless ';
    for k=1:length(mach_optlist),
      msg = [msg, sprintf('''%s'', ',mach_optlist{k})];
    end
    msg=[msg, 'option(s) on Tecmag'];
    warning(msg);
  end
  blnktail=4; %% trailing 4 channels are 0 on TMG
else %%Stelar, Parma
  blnktail=0; %% no trailing blank on Stelar and others
end

if ~(t_fl | f_fl), 
  f_fl=1; t_fl=1; %%default condition: T_domain & FFT
end 
  
if t_fl, %%if time domain
  fname=[NMRpath,fname];
  if (I_fl),
    TDfids_o=TDfids; %% incremental loading : backup fids and seqs 
    FDfids_o=FDfids;
    NMRseq_o=NMRseq;
    NMRamp_o=NMRamp;
    NMRshift_o=NMRshift;
    NMRnoise_o=NMRnoise;
    NMRmisc_o=NMRmisc;
  end

  if (Instrument == 1), %%Stelar
    mach_optlist=cat(2,{'-m'},mach_optlist); %%expand multiblock to multizone  
    [TDfids,pseq,NMRmisc]=slrload(fname, mach_optlist{:});
    pseq(5,:)=[]; %%multiblock expanded to multizone, NBLK was 1 
  elseif (Instrument == 2) %%Parma
    mach_optlist=cat(2,{'-n','-w'},mach_optlist); %%Normalize, shift workaround
    [TDfids,pseq,NMRmisc]=pnaload(fname, mach_optlist{:});
    pseq(5,:)=[]; %%remove 
  elseif (Instrument == 3) %%Ac300
    mach_optlist=cat(2,{'-n','-y'},mach_optlist); %%Normalize, shift Ydata
    [TDfids,pseq,NMRmisc]=ac300load(fname, mach_optlist{:});
    NMRmisc{1}=cellify(NMRmisc{1});
  else  %%TecMag
    [TDfids,pseq,NMRmisc]=tmgload(fname,'-n','-t'); 
    NMRmisc=cellify(NMRmisc);
  end
  sz=size(TDfids); 
  szps=size(pseq);
  if (szps(1)< 9),
    S_fl = 0; %%less than 3 pulses, not SSE
  end

  %%NMRseq=[N_td;dwell_td;freq;Nscans; 
  %%acqd;chan_topecho;N_fd;dwell_fd,delay1,..]
  NMRseq=zeros(8 + floor((szps(1)-4)/2),sz(2));
 
  NMRseq(1:4,:)=pseq(1:4,:);
  if (szps(1)>4), 
    NMRseq(5,:)=pseq(szps(1),:);
  end
    
  if (szps(1)>6),
    if (S_fl), 
      NMRseq(6,:)=(.5*pseq(szps(1)-5,:) + pseq(szps(1)-4,:) +...
	  .5*(pseq(szps(1)-3,:)-pseq(szps(1)-1,:)) - ...
	  pseq(szps(1),:))./pseq(2,:); 
    else
      NMRseq(6,:)=(.5*pseq(szps(1)-3,:) + pseq(szps(1)-2,:) - ...
	  pseq(szps(1),:))./pseq(2,:); 
    end
  end

  for k=5:2:szps(1);
    NMRseq(9 + (k-5)/2,:) = (pseq(k,:) > 0) .* (.5*pseq(k,:) + pseq(k+1,:));
    if ((k+2) <= szps(1)),
      NMRseq(9 + (k-5)/2,:) = NMRseq(9 + (k-5)/2,:) + ...
	  (pseq(k+2,:) > 0) .* pseq(k+2,:) /2;
    end
  end
end %%END if t_fl
%%

sz=size(TDfids);
if (sz(2)>1),
  ix_=sum((NMRseq([1:2 4:6],2:end)-NMRseq([1:2 4:6],1:end-1)) ~= 0);
  hx=[0 find(ix_) sz(2)]; %%hx delimits blocks of homogeneous zones
else
  hx=[0 1];
end

if (f_fl), 
  NMRamp=zeros(1,sz(2));
  NMRshift=zeros(1,sz(2));  
  %%set FFT base, enlarge fft buffer if zero padding is enabled
  fsz=sz;

  if (n_fl),
    NMRnoise=zeros(1,sz(2));
  else
    NMRnoise=1e3./sqrt(NMRseq(4,:));
  end
  
  if (z_fl),  
    nn_ = 2^ceil(log(abs(z_fl))/log(2)); 
    if (z_fl < 0 | nn_ <= sz(1)),
      FDfids=zeros(sz(1),sz(2)); 
      NMRseq(7,:)=sz(1)*ones(1,sz(2));
    else %%zero paddind up to nn_ points
      FDfids=zeros(nn_,sz(2)); 
      fsz(1)=nn_; %%pad fid buffer with zeroes
      NMRseq(7,:)=nn_*ones(1,sz(2));
    end
  else
    FDfids=zeros(sz(1),sz(2)); 
    NMRseq(7,:)=NMRseq(1,:);
  end  
  NMRseq(8,:)=1../(NMRseq(7,:).*NMRseq(2,:)); %%vide infra
  
  if (T_fl(1)), %% interactive apod. window
    gcf; %%xselect();
    [curs,apod_p]=xfidfilt('ria',1);
    k_fl=1; pblank=curs(1);
    T_fl(2)=1; Twin=curs(2);
    b_fl=1; b_lim=[curs(3) inf];
  end

  %%Preparing for FFT
  pblk_=1; 
  win__bck=zeros(length(hx)-1,2);
  blc__bck=zeros(1,sz(2));

  for ii=1:length(hx)-1; %%LOOP over EXPERIMENTS 
    ixs=hx(ii)+1:hx(ii+1); ix=ixs(1);
    N=NMRseq(1,ix);	%%n. of meaniningful channels 
    %%NMRseq(5,ix)=acq. delay
  
    if (~T_fl(2)),
      Twin=NMRseq(6,ix); fl_=1; %% use chan_topecho; fl_==0 -> use time units
    else
      fl_=T_fl(2)-1;
    end
    win_=time2ch(Twin,ix,fl_);  
    if (length(win_) < 2),
      win_=[win_,1../(apod_p.fwid*NMRseq(2,ix))]; %%win width: use default
    end
    win__bck(ii,:)=win_;
    
    if (b_fl),   %%if baseline correction
      blim_=time2ch(b_lim,ix,b_fl-1); %%conversion to channel units
      blim_=round(sort(blim_));
      blim_(1)=max(1,blim_(1)); %%%bug in Scilab code
      blim_(2)=min(N-blnktail,blim_(2));
      if ((blim_(2)-blim_(1)) >= 1),
        blc = mean(TDfids(blim_(1):blim_(2),ixs)); 
      elseif ((blim_(2)-blim_(1)) < 0),
        blc = zeros(size(ixs)); 
      else
        blc = TDfids(blim_(1),ixs); 
      end
    else 
      blc=zeros(size(ixs));
    end 
    blc__bck(ixs)=blc; %% store baselines for future usage
    
    if (k_fl),   %%echo preblanking: blank 1:pblk_-1 channels
      pblk_=ceil(time2ch(pblank,ix,k_fl-1)); %%conversion to channel units
    end 

    %% apodization window
    apw=eval(sprintf(...
      '%s(N-blnktail,win_(1),win_(2),apod_p.half)',apod_p.typ))';
    L=NMRseq(7,ix);
    fidb=zeros(L,length(ixs)); 
    rng=(1:N-blnktail)'; %%TD range: Tecmag yields last 4 chans = 0
    srng=fix(1+mod(rng-win_(1),L));	%shifted range: top of echo on 1st
                                     %channel %%bug in Scilab code				
    %% apodize and do FFT
    fidb(srng,1:end) = (TDfids(rng,ixs) - ones(size(rng))*blc) .* ...
	(((rng >= pblk_).*apw) * ones(size(ixs)));
    f_fidb=fft(fidb);
    rng=1+rem(L/2+(0:L-1),L)'; 

    FDfids(1:L,ixs)=f_fidb(rng,:); %%fftshift

  end %%End loop over experiments

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% FFT phasing (if enabled)
  RealFFT=0;
  if ((length(Phcf) == 1) & Phcf),
    Phcf=xphase(1,'-n');
    RealFFT=1;
    if (min(size(Phcf))==1) % homo-phasing
      Phcf=Phcf(:)'; %force row
    elseif (all(all(~diff(Phcf')))), %phase coeff are const over zones
      Phcf=Phcf(:,1)';
    end
  elseif (length(Phcf) > 1),
    RealFFT=1;
    if (size(Phcf,2) > 1),
      if (length(Phcf) ~= 2  & size(Phcf,2) ~= size(NMRseq,2)),
        disp('Warning: previous phasing cannot be applied. Size mismatch');
        %phasing failure, but use real(FFT) anyway
	%%
	%%keyboard;
	%%
      else
        FDfids=rephase(FDfids,Phcf,NMRseq(7,:));
      end
    end
  end
        
    
  if (F_fl(1)), %% interactive seletion of FFT window
    gcf; %%xselect();
    Fwin=xfftsel([0 0 1]+RealFFT*[1 0 -1],1);
    F_fl=[0 1];
  end
  

  %% evaluation of amplitudes and noise
  for ii=1:length(hx)-1; %%RE-LOOP over EXPERIMENTS 
    ixs=hx(ii)+1:hx(ii+1); ix=ixs(1);
    N=NMRseq(1,ix);  %%n. of meaniningful TD channels 
    L=NMRseq(7,ix);  %%n. of meaniningful FD channels 

    fwin_=freq2ch(Fwin,ix,F_fl(2)-1);
    if (length(fwin_)<2), %% use default fft width
      fwin_=[fwin_,freq2ch(NMRfftwin('fwid'),ix,NMRfftwin('funit'))];
    end

    if (fwin_(1)), %% shift of FFT: this affects FDfids 
      %%dephasing 
      dp=(fwin_(1)*2*pi*[0:(L/2-1), -(L/2):-1]/L)'*ones(size(ixs)); 
      %%frequency shift
      FDfids(:,ixs)=ifft(fft(FDfids(:,ixs)).*exp(i*dp)); 
    end

    ft_rng=max(1,1+ceil(L/2-fwin_(2))):min(L,1+floor(L/2+fwin_(2)));
    if (length(ft_rng) > 1),
      if (RealFFT),
        [t__,iy]=max(abs(real(FDfids(ft_rng,ixs)))); 	%%iy:relative index
        for jj=1:length(ixs), %take real(FFT) at max(abs(real(FFT))))           
          t__(jj) = real(FDfids(iy(jj)+ft_rng(1)-1,ixs(jj))); 
        end  
      else
        [t__,iy]=max(abs(FDfids(ft_rng,ixs))); 	%%iy:relative index
      end
    else
      if (RealFFT),
        t__=real(FDfids(ft_rng,ixs));
      else
        t__=abs(FDfids(ft_rng,ixs));
      end
      iy=ones(size(ixs));
    end
    iy=iy+ft_rng(1)-1;  %% index of t__ in FDfids(:,ix)
    NMRamp(ixs)=t__;   %% maximum within fft window    
    NMRshift(ixs)=(mod(fwin_(1)+iy-1,L)-L/2)*NMRseq(8,ix);

    if (n_fl), %%estimate of noise from the noisy echo tail 
      if (b_fl),
        %%use baseline corr. window as noise win
        nn_=max(2,time2ch(b_lim(1),ix,b_fl-1)); 
      else
        nn_=max(2,round(2/3*L)); %%default noise win: [2/3 1]
      end
      fidb=zeros(L,length(ixs)); 
      n_rng=(1:(N-nn_-blnktail+1))';

      if (~length(n_rng)),
	NMRnoise(ixs)=1e3./sqrt(NMRseq(4,ixs)); %conventional value
      else
	apw=eval(sprintf(...
	    '%s(N-blnktail,win__bck(ii,1),win__bck(ii,2),apod_p.half)',...
	    apod_p.typ))';
	nn_=fix(nn_);
	%%
	fidb(n_rng,1:end)=(TDfids(n_rng+nn_-1,ixs) - ...
			   ones(size(n_rng))*blc__bck(ixs)).*...
	    (apw(n_rng)*ones(size(ixs))) * sum(apw)./sum(apw(n_rng));
	j0 = (fwin_(1)+iy - L/2) == 1;
	ix0 = ixs(j0); %%detection at D_nu==0
	j1 = (fwin_(1)+iy - L/2) ~= 1;
	ix1 = ixs(j1); %%detection at D_nu~=0
	if (length(n_rng)>1),
	  NMRnoise(ix0) = sqrt(sum(conj(fidb(n_rng,j0)).*fidb(n_rng,j0)));
	elseif (length(n_rng)==1),
	  NMRnoise(ix0) = sqrt(conj(fidb(n_rng,j0)).*fidb(n_rng,j0));
	end
	if (length(ix1)),
	  f_fidb=fft(fidb(:,j1));
	  dp=(fwin_(1)*2*pi*(rng-L/2-1)/L) * ones(size(ix1)); 
	  f_fidb=ifft(fft(f_fidb).*exp(i*dp)); %%frequency shift
	  %%trick: NMRnoise (matrix) is accessed like a vector
	  NMRnoise(ix1)=abs(f_fidb(1+mod(iy(j1)-1-L/2,L) + ...
				   (0:length(ix1)-1)*L))'; %%noise
	end
      end
    end %%end if n_fl

  end %%End re-loop over experiments
end %% END if f_fl
  
if (I_fl & t_fl), %%Incremental loading: join old and new data
  sz__=size(TDfids_o);
  d_trw=sz__(1)-sz(1);
  d_frw=size(FDfids_o,1)-fsz(1);
  d_srw=size(NMRseq_o,1)-size(NMRseq,1);

  TDfids=[[TDfids_o;zeros(-d_trw,sz__(2))],[TDfids;zeros(d_trw,sz(2))]]; 
  FDfids=[[FDfids_o;zeros(-d_frw,sz__(2))],[FDfids;zeros(d_frw,sz(2))]]; 
  NMRseq=[[NMRseq_o;-1*ones(-d_srw,sz__(2))],[NMRseq;-1*ones(d_srw,sz(2))]]; 
  NMRamp=[NMRamp_o,NMRamp];
  NMRshift=[NMRshift_o,NMRshift];
  if (size(NMRnoise_o) == size(NMRamp_o)) 
    NMRnoise=[NMRnoise_o,NMRnoise];
  else
    disp(['Errobars of previous dataset have been manually ',...
	  'altered, you cannot append errorbars to NMRnoise']);
  end

  %% merge NMRmisc row-wise for tnt and column-wise for PNA, SLR
  wrong=0;
  if (Instrument==1 | Instrument==2), %%Stelar, PNA
    if (iscell(NMRmisc_o) | (size(NMRmisc,1) ~= size(NMRmisc_o,1))),
      wrong=1;
    else
      NMRmisc=cat(2, NMRmisc_o, NMRmisc);
    end
  elseif (~Instrument), %Tecmag
    if (~iscell(NMRmisc_o) | iscell(NMRmisc_o{1})),
      wrong=1;
    else
      NMRmisc=cat(1, NMRmisc_o, NMRmisc); 
    end
  else %AC300
    if (~iscell(NMRmisc_o) | ~iscell(NMRmisc_o{1}) | ...
        (size(NMRmisc)~=size(NMRmisc_o)) | ...
        (length(NMRmisc)==2 & size(NMRmisc{2},1)~=size(NMRmisc_o{2},1))),
      wrong=1;
    else
      NMRmisc{1}=cat(1, NMRmisc_o{1}, NMRmisc{1}); 
      if (length(NMRmisc)==2),
        NMRmisc{2}=cat(2, NMRmisc_o{2}, NMRmisc{2}); 
      end
    end
  end
  if (wrong),
    disp('Miscellaneous data NMRmisc are heterogeneous, cannot be merged');
  end
end 

%%find out which is the relevant 2D par (default is frequency)
NMRpar2d = struct('name2d','','val',[]);
[NMRpar2d.val,NMRpar2d.name2d]=tell2d(0);

%%save status
T_fl(1)=0;
NMRstat={b_fl,b_lim,k_fl,pblank,n_fl,z_fl, T_fl,Twin,F_fl,Fwin, apod_p, Phcf}; 
NMRnoise=real(NMRnoise);


%%PRIVATE functions

function r=isblank(s);
r=(~length(s) | all(s==' '));

function [v,chfl]=getopt(str,n);
v=[];
if (length(str)<3),
  chfl=0;
  return;
else
  chfl=(str(3)=='c');
  if (chfl & length(str)<4),
    return;
  end
end
ss=str((3+chfl):length(str));
if isblank(ss),
  return
end
v=sscanf(ss,'%g:').'; 
v=v(1:min(length(v),n));


function ch=time2ch(t,k,chfl)
global NMRseq;
if (chfl), 
  ch=t; 
else 
  ch=t/NMRseq(2,k); 
end

function ch=freq2ch(f,k,chfl);
global NMRseq;
if (chfl), 
  ch=f; 
else 
  ch=f/NMRseq(8,k); 
end

%% Copyright (C) G. Allodi
%% Revision 2.1   20-may-2004
