function nmrdisplay(varargin); 
%% interactive display of time-domain or frequency domain fids
%% options: 
%%  -f: frequency domain (default: time domain)
%%  -m[ria] : real, imag, abs (default:ria)
%%  -n: channel units istead of sec/Hz
global TDfids FDfids NMRseq;

%%defaults:
ria_fl=[1 1 1];
tsca=1;
FD=0; %% time domain

for k=1:nargin;
  ar=varargin{k}; l_ar=length(ar);
  if (~isstr(ar)), 
    error('All input arguments must be strings');
  end
  swtc = ar(1:min(2,l_ar));
  opt = ar(3:l_ar);
  %% Check options
  if (strcmp(swtc,'-m')),
    ria_fl(1) = (0 < length(find(opt=='r')));
    ria_fl(2) = (0 < length(find(opt=='i')));
    ria_fl(3) = (0 < length(find(opt=='a')));
    if (~any(ria_fl)),
      disp(sprintf(['Valid options are -mr, -mri, -mria, etc.\n ',...
	  'Applying default.\n']));
      ria_fl=[1 1 1];
    end
  elseif (strcmp(swtc,'-n')),
    tsca=0;
  elseif (strcmp(swtc,'-f')),
    FD=1;
  else
    warning(sprintf(' unknown ''%s'' option.',swtc));
  end
end
  
%% program setup
[par2d,name2d,varseqFL]=tell2d(FD);
ufacts=[1e6 1e-6]; %% unit conversion factors: [s Hz] -> [us, MHz]

cix=1;
do_init=1;
funcs={'real','imag','abs'};
clrs=['r','g','b'];

figure(gcf);
clf;


%%button names 
guiname={'Curs#1','Curs#2','Blowup','Unzoom','View','+8 zone','+1 zone',...
	 '-1 zone','-8 zone','Done'};
zinc=[8,1,-1,-8]; %% zone incremnts

pos=[.08,.0,.08,.06]; 
dxpos=.14;
dypos=.115;
dpos=[dxpos,0];

% instal GUIs
for ii=1:prod(size(guiname));
  clbck=sprintf('set(gcf,''userdata'',%g);uiresume(gcf)',ii);  
  if (ii==3),
    dpos=[0,-dypos];
    pos(1:2)=[.915,pos(2)+(length(guiname)-3)*dypos];
  end
  ui=uicontrol('style','push','units','normal','pos',pos,'string', ...
    guiname{ii},'call',clbck,'Tag','UINMR');
  pos(1:2)= pos(1:2)+dpos;
end

%loop until Done button pressed (async polling)  
while (1),
  %%main loop
  if (do_init),
    if (FD),
      mtxname='FDfids';
    else
      mtxname='TDfids';
    end
 
    %%initialize variables for 1st plot
    nn=1:NMRseq(1+6*FD,cix);
    nnlims=[nn(1) nn(end)]; %% Used by 'blow up'
    xx=nn - FD*(nn(end)/2 + 1);  %%shift for FD: FFT is  centred around 0  
    if (tsca),
      xx = xx*NMRseq(2+6*FD,cix)*ufacts(1+(FD~=0));
    end
    xxlims=[xx(1) xx(end)];
    crs=xxlims;
    xpos=crs'*[1 1];
    do_init=0;

    %first plot		
    pid=-1*ones(1,3);
    max_y=0;
    flg=1;
    for ii=3:-1:1;
      if (ria_fl(ii)),
	yy=eval(sprintf('%s(%s(nn,cix))',funcs{ii},mtxname));
	if (flg),
	  if (ii==3), %abs
	    max_y=max(yy);
	    flg=0;
	  else
	    max_y=max(max_y,max(abs(yy)));
	  end
	end
	pid(ii)=plot(xx,yy,clrs(ii));
	hold on;
      end
    end
    epid=zeros(1,3); %handles of baseline & cursors: initialize
    epid(3)=line([xx(1) xx(end)],zeros(1,2),'color','k','linestyle','--');
    if (tsca),
      if (FD),
	xlabel('Frequency shift (MHz)');
      else
	xlabel('time (\mus)');
      end
    else
      xlabel('channels');
    end
    ylabel('amplitude (a.u.)');
    ax=axis; ax(1:2)=[xx(1) xx(end)]; axis(ax);
    epid(1)=line(xpos(1,:),ax(3:4),'color','b');
    epid(2)=line(xpos(2,:),ax(3:4),'color','m');
    hold off
    titid=title(sprintf('Zone %d, %s = %g',cix,name2d,par2d(cix)),...
		'fontsize',12);
    set(gcf,'userdata',0); 
    blowup=0;
  end  %%if do_init

  uiwait(gcf);  
  %callback definitions for gui buttons
  jj=get(gcf,'userdata');

  if (jj < 3), %% moving a cursor 
    set(epid(jj),'linestyle','--');
    mpos=ginput(1);
    crs(jj)=min(max(xx(1),mpos(1)),xx(end));
    xpos(jj,:)=crs(jj)*[1 1];    
    set(epid(jj),'linestyle','-','xdata',xpos(jj,:));
  elseif (jj==3),  %%blow up
    crs=sort(crs); 
    [nnlims xxlims]=crs2lims__(crs,tsca,cix,FD); %%recalc limits
    blowup=1;
    z_=axis; z_(1:2)=xxlims; axis(z_);
  elseif (jj==4), %%reset blow up
    if (blowup),
      %%default limits
      [nnlims xxlims]=crs2lims__([xx(1) xx(end)],tsca,cix,FD); 
      z_=axis; z_(1:2)=xxlims; axis(z_);      
      blowup=0;
    end
  elseif (jj==5)
    [ria_,sca_,fd_]=plotmodesel_(ria_fl,tsca,FD);
    if (length(fd_)),
      ria_fl=ria_;
      FD=fd_;
      tsca=sca_;
      do_init=1;
    end
  elseif (jj < 10), %%incrementing zone
    cix=max(min(cix+zinc(jj-5),size(TDfids,2)),1);
    if (varseqFL)
      nn=1:NMRseq(1+6*FD,cix);
      xx=nn;
      if (FD),
	xx=xx - nn(end)/2 - 1;
      end
      if (blowup), %%validate limits
	nnlims=min([nnlims,NMRseq(1+6*FD,cix)*[1 1]]);
      else
	nnlims=[1 nn(end)];
      end
      xxlims=nnlims-FD*(nn(end)/2 + 1);
      if (tsca),
	xx=xx*NMRseq(2+6*FD,cix)*ufacts(1+(FD~=0));
	xlims=xlims*NMRseq(2+6*FD,cix)*ufacts(1+(FD~=0));
      end
      for k=1:2;
	crs(k)=min(max(xx(1),crs(k)),xx(end)); %validate cursors if out
                                               %of bounds
        xpos(k,:)=crs(k)*[1 1];                
        set(epid(k),'xdata',xpos(k,:),'ydata',ax(3:4));
      end
      [ax,max_y]=replottf(pid,epid,xx,cix,FD);
      set(epid(2),'xdata',[xx(1) xx(end)]);
      %%validate cursors (if out of bounds)
      for k=1:2;
	xpos(k,:)=ascis(crs(k),cix); crs(k)=xpos(k,1);
      end
    else
      [ax,max_y]=replottf(pid,epid,[],cix,FD);
      for k=1:2;
        set(epid(k),'ydata',ax(3:4));
      end      
    end
    set(epid(1:2),'visible','on');
    titid=title(sprintf('Zone %d, %s = %g',cix,name2d,par2d(cix)),...
		'fontsize',12);
  elseif (jj==10), %done
    break;
  end
end

%clean-up: delete buttons
delete(findobj('Tag','UINMR'));


function [nnl,xxl]=crs2lims__(crs,sca,cix,fd)
%%new limits and cursors for blowup/reset
global NMRseq;

ufacts=[1e-6 1e-6];
ceps=eps*max(100,NMRseq(1+6*fd,cix));
if (sca),
  nnl=crs/(NMRseq(2+6*fd,cix)*ufacts(1+(fd~=0)));
else
  nnl=crs;
end
%%round curs position and shift if in FFT mode
nnl=[floor(nnl(1)+ceps) ceil(nnl(2)-ceps)];
xxl=nnl;
if (fd),
  nnl=nnl + NMRseq(1+6*fd,cix)/2 + 1;
end
%%nnlims=[max(1,floor(nnlims(1))) min(nn(end),ceil(nnlims(2)))]
if (sca),
  xxl=xxl*NMRseq(2+6*fd,cix)*ufacts(1+(fd~=0));
end


%%% Extern:
%%% function [ria_,sca_,fd_]=plotmodesel_(ria,sca,fd);

%(private) REPLOTTF: replots fid from a new zone
function [ax,max_y]=replottf(pid,epid,xvec,zone,fd);
%pid=[id_real,id_imag,id_abs];
%epid=[handles of cursors, etc.];
%if xvec==[], xdata kept 
%zone: new zone index (no check whether it exists)
%fd~=0: frequency domain;
global NMRseq;
if (fd),
  global FDfids;
  mtxnam='FDfids';
else
  global TDfids;
  mtxnam='TDfids';
end

set(epid(1:2),'visible','off'); %I don't want these to affect scaling
N=NMRseq(1+6*fd,zone);
fcn={'real','imag','abs'};
max_y=0;
fl=1;
%seek max_abs
for ii=3:-1:1;
  if (pid(ii)>0)
    y=eval(sprintf('%s(%s(1:N,zone))',fcn{ii},mtxnam));
    if (length(xvec)),
      set(pid(ii),'Xdata',xvec,'Ydata',y);
    else
      set(pid(ii),'Ydata',y);
    end
    if (fl),
      if (ii==3),
        max_y=max(y);
        fl=0;
      else
        max_y=max(max_y,max(abs(y)));
      end
    end
  end
end
if (~length(xvec)),
  xvec=get(pid(ii),'Xdata');
end

axis auto;
ax=axis; ax(1:2)=[xvec(1),xvec(end)]; axis(ax);




