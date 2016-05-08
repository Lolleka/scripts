function phcf=xphase(varargin); 
%% arguments:
%%   phcf_init = [phase1, phase2]: initial phase coefficients 
%%   fsca (flag): if non zero, use time scale in Hz instead of channels
%%   option switches (strings): 
%%           '-f': call xfftsel at exit and re-evaluate NMRamp from real(FTfids);
%%           '-n': do not re-evaluate NMRamp;
%%           '-q': query: re-evaluate NMRamp?
%%           '-m': mark all zones;
%% output:
%%   phcf -> phasing coefficients (for all zones, except those marked) 

global FDfids NMRseq;

fsca=0; %% defaults
phcf=[0; 0]; 
exitmode=1; %exitmode=1 -> '-q' option; exitmode=0 -> '-n'; exitmode= 2 -> '-f'

%%phasing coeff buffer 
cf_buff=zeros(3,size(NMRseq,2)); %% [[ph1,ph2;flag],...]
status=['      ';'MARKED'];

%% parsing input arguments 
for ii=1:nargin;
  if (ispair(varargin{ii})),
    phcf=varargin{ii}(:);
  elseif (isscal(varargin{ii})),
    fsca=varargin{ii};
  elseif (isstr(varargin{ii})),
    switch varargin{ii},
    case '-f', exitmode=2;
    case '-q', exitmode=1;
    case '-n', exitmode=0;
    case '-m', cf_buff(3,:)=ones(1,size(NMRseq,2)); %mark all zones
    otherwise 
      warning(sprintf('%s: invalid option switch "%s"',mfilename,varargin{ii}));
    end
  else
    warning(sprintf('%s: invalid argument #%d',mfilename,ii));
  end
end
phcf_init=phcf;

%% program setup
[par2d,name2d,varseqFL]=tell2d(1);
zinc=[8,1,-1,-8]; %% zone incremnts
if isunix,
  fntsz=12;
  uifntsz=12;
else
  fntsz=10;
  uifntsz=9;  
end
%%button/slider  names 
guiname={'+8 Zone','+1 Zone','-1 Zone',...
	 '-8 Zone','Mark','Reset','Done','Zoom+','Zoom-','Unzoom'};

%%initialize plot
zoomx=0;
cix=1;
nn=1:NMRseq(7,cix);
fd_fid=rephase(FDfids(1:nn(end),cix),phcf,nn(end));

if (fsca),  
  xx = (nn - nn(end)/2 - 1)*NMRseq(8,cix); 
  d_x=NMRseq(8,cix);
  x0=d_x*(nn(end)/2+1);
else
  xx = nn - nn(end)/2 - 1;
  x0 = nn(end)/2 + 1;
  d_x=1;
end

%first plot		
clrs=['r','b'];
%funcs={'real','imag'};
epid=zeros(1,2); %handles of cursors: initialize
pid=plot(xx,real(FDfids(nn,cix)),clrs(1),xx,imag(FDfids(nn,cix)),clrs(2));
max_y=max(abs(FDfids(nn,cix)));
if (fsca),
  xlabel('Frequency (Hz)');
else
  xlabel('channels');
end
ylabel('amplitude (a.u.)');
axbox=[xx(1) xx(end) -max_y max_y]; axis(axbox); %force axes to match figure 
epid(1)=line(axbox(1:2),zeros(1,2),'color','k','linestyle','-');
epid(2)=line(zeros(1,2),axbox(3:4),'color','k','linestyle','-');

titid=title(sprintf('Zone %d, %s = %g  %s',...
  cix,name2d,par2d(cix),status(cf_buff(3,cix)+1,:)),'fontsize',fntsz);
udt=struct('ui',0,'val',phcf,'txt',zeros(1,3),'sli',[0 0]); 

%% install GUIs
%1)buttons
pos=[.02,.0,.08,.056]; 
dxpos=.095;
dpos=[dxpos,0];
for ii=1:prod(size(guiname));
  clbck=sprintf(['u=get(gcf,''userdata'');u.ui=%d;set(gcf,''userdata'',u);',...
		 'uiresume(gcf)'],ii);
  ui=uicontrol('style','push','units','normal','pos',pos,'string', ...
	       guiname{ii},'call',clbck,'Tag','UINMR');
  if (ii==4),
    pos(1)=.58;
  elseif (ii==7),
    pos(1:2)=[.02-dxpos,1-pos(4)];
  end
  pos(1:2)= pos(1:2)+dpos;
end

%2)sliders & text
pos=[.915 .055 .035 .88];
tpos=[.135 .875 .12 .045];
lms=[pi pi/2^max(3,ceil((log(nn(end)/512))/log(2)))];
for ii=((1:2)+length(guiname));
  jj=ii-length(guiname);
  clbck=sprintf(['u=get(gcf,''userdata'');u.ui=%d;v=get(gcbo,''value'');' ...
		 'u.val(%d)=v;set(gcf,''userdata'',u);set(u.txt(%d),' ...
		  '''string'',num2str(v*180/pi));uiresume(gcf)'],ii,jj,jj);
  
  udt.sli(jj)=uicontrol('Style','slider','Units','normalized','Position',...
			pos,'Value',0,'Max',lms(jj),'Min',-lms(jj),...
			'sliderstep',[.001 .025]*lms(jj),'callback',clbck, ...
			'Tag','UINMR');
  
  u=uicontrol('Style','text','Units','normalized','Position',...
	      pos+[-.013 .89 .013 -.84] ,'Tag','UINMR','string',...
	      sprintf('ph%d',jj),'fontsize',uifntsz);
  pos(1)=pos(1)+.05;
  %text
  udt.txt(jj)=uicontrol('Style','text','Units','normalized','Position',...
			tpos,'Tag','UINMR','string',num2str(phcf(jj)*...
			180/pi),'fontsize',uifntsz);
  tpos(1)=tpos(1)+.645;
end
%text: zooming level
udt.txt(3)=uicontrol('Style','text','Units','normalized','Position',...
  [.02,.88,.065,.045],'Tag','UINMR','string','X 1',...
  'fontsize',uifntsz);

set(gcf,'userdata',udt);


figure(gcf);
phcf_bck=phcf;
%loop until Done button pressed (async polling)  
while(1),
  uiwait(gcf);  
  %callback definitions for gui buttons
  udt=get(gcf,'userdata');
  replotfl=1;
  syncfl=0;
  if (udt.ui <5), %% incrementing zone
    bckix=cix;
    cix=max(min(cix+zinc(udt.ui),size(FDfids,2)),1);
    nn=1:NMRseq(7,cix);
    if (varseqFL)  %%sequence varies throughout zones: re-initializing
      if (fsca),  
        xx = (nn - nn(end)/2 - 1)*NMRseq(8,cix); 
      else
        xx = nn - nn(end)/2 - 1;
      end
    end
    syncfl=1;
    if (cf_buff(3,cix)), %%if marked: restore previous phasing
      phcf=cf_buff(1:2,cix);
    else
      phcf=phcf_bck;
    end
  elseif (udt.ui==5),  %Mark/Unmark
    replotfl=0;
    cf_buff(3,cix)=~cf_buff(3,cix); %%toggle flag
    if (cf_buff(3,cix)),
      cf_buff(1:2,cix)=phcf(:); %%store
    else
      phcf_bck=phcf; %%on unmarking, make current phase the default
    end
  elseif (udt.ui==6),  %Reset
    syncfl=1;
    resp=questdlg('Reset phase of:','xphase dialog','All zones',...
		  'Unmarked only','Cancel','Cancel');
    if (~strncmp(resp,'Can',3)),
      phcf_bck=phcf_init(:);
      if (strncmp(resp,'All',3)),
	phcf=phcf_bck;
	cf_buff=zeros(size(cf_buff));
      elseif (~cf_buff(3,cix)), %if 'Unmarked'
	phcf=phcf_bck;
      end
    end
  elseif (udt.ui==7), %%Done
    resp=questdlg('Do you really want to exit?','xphase dialog','Yes',...
      'No','No');
    if (strcmp('Yes',resp)),
      break;
    end
  elseif (udt.ui<=9), %zoom+-
    zoomx=zoomx-sign(udt.ui-8.5);
    zoomx=max(0,min(round(log(NMRseq(7,cix))/log(2)-2),zoomx));
    axis([axbox(1:2)/2^zoomx, axbox(3:4)]);
    set(udt.txt(3),'string',sprintf('X %d',2^zoomx));
  elseif (udt.ui==10), %unzoom
    zoomx=0;
    axis(axbox);
    set(udt.txt(3),'string','X 1');
  elseif (udt.ui>length(guiname)), %sliders
    k=udt.ui-length(guiname);
    v=get(udt.sli(k),'value');
    phcf(k)=v;
    if (cf_buff(3,cix)), %%if marked
      cf_buff(1:2,cix)=phcf(:); %%store
    else
      phcf_bck=phcf; %update default phase setting
    end    
  end 
  
  if (syncfl), %sync sliders to phcf
    for ii=1:2;
      set(udt.sli(ii),'value',phcf(ii));
      set(udt.txt(ii),'string',num2str(phcf(ii)*180/pi));
    end
  end
  if (replotfl), %rephase and replot
    fd_fid=rephase(FDfids(1:nn(end),cix),phcf,nn(end));
    if (varseqFL & (udt.ui<5)), %zone was incremented: redraw axes
      set(pid(1),'xdata',xx,'ydata',real(fd_fid));
      set(pid(2),'xdata',xx,'ydata',imag(fd_fid));
      set(epid(1),'xdata',[xx(1) xx(end)]);
    else
      set(pid(1),'ydata',real(fd_fid));
      set(pid(2),'ydata',imag(fd_fid));
    end
    if (udt.ui<5),
      max_y=max(abs(FDfids(:,cix)));
      axbox=[xx(1) xx(end) -max_y max_y];
      set(epid(2),'ydata',axbox(3:4));
      axis([axbox(1:2)/2^zoomx, axbox(3:4)]);
    end
  end
  if (udt.ui<6),
    titid=title(sprintf('Zone %d, %s = %g  %s',cix,name2d,...
      par2d(cix),status(cf_buff(3,cix)+1,:)),'fontsize',fntsz);
  end
end %end while

%clean-up: delete buttons
delete(findobj('Tag','UINMR'));

if (any(phcf) | any(sum(abs(cf_buff(1:2,:))).*cf_buff(3,:))),
  resp=questdlg({'Quitting','Which zones are to be rephased?',...
      '(NB: FDfids is affected)'},'xphase dialog','All zones',...
    'Marked only','None','All zones');
  atexit=find(strncmp({'Al','Ma','No'},resp,2));
else
  atexit=3;
end

if (atexit < 3), %%exit without 'quit'
  mtx_flg=0;
  if (any(cf_buff(3,:))), %%there are marked zones, output a matrix 
    phcf=zeros(2,size(FDfids,2));
    mtx_flg=1;
  end
  for k=1:size(FDfids,2);
    if (cf_buff(3,k)), %%marked zone
      FDfids(:,k)=rephase(FDfids(:,k),cf_buff(1:2,k),NMRseq(7,k));
      phcf(:,k)=cf_buff(1:2,k);
    elseif (atexit==1), %%unmarked, apply default phasing
      FDfids(:,k)=rephase(FDfids(:,k),phcf_bck,NMRseq(7,k));
      if (mtx_flg),
        phcf(:,k)=phcf_bck;
      end        
    end
  end  
else
  phcf=[0; 0];
end

if (exitmode==1), %%query: re-evaluate NMRamp
  resp=questdlg('Evaluate amplitude from rephased real(FFT(fid))?',...
    'xphase dialog','Yes','No','Yes');
  exitmode=strcmp(resp,'Yes');
end 

if (exitmode),
  NMRamp=zeros(1,size(NMRseq,2));
  crs=xfftsel([1 0 0],0,1);
  crs=[crs(1)-crs(2), crs(1)+crs(2)];
  for ii=1:size(NMRseq,2);
    rng=round([max(-NMRseq(7,ii)/2,crs(1)), min(NMRseq(7,ii)/2,crs(2))] + NMRseq(7,ii)/2);
    [dmy,jx]=max(abs(real(FDfids(rng(1):rng(2),ii))));
    NMRamp(ii)=real(FDfids(rng(1)+jx-1,ii));
  end
end
  
if (~exitmode & (atexit==3)),
  disp([mfilename,': nothing done'])
end


%%local macros
function r=isscal(a);
r=isnumeric(a) & ~isempty(a) & (max(size(a))==1);

function r=ispair(a);
r=isnumeric(a) & (length(a)==2);

