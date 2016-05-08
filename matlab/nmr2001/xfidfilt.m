function [crs,apod_p]=xfidfilt(ria_fl,tsca); 
% arguments:
%	ria_fl=[real_flag,imag_flag,abs_flag] 
%      tsca=1 -> uses time scale il secs instead of channels
global TDfids NMRseq NMRapod;

% validation of arguments, apply defaults
if (nargin < 1),
  ria_fl=[1 1 1];
elseif (size(ria_fl(:),1) < 3),
  ria_fl=[1 1 1];
elseif (~any(ria_fl))
  ria_fl=[1 1 1];
end;
if (nargin < 2),
  tsca=0;
end

% program setup
apod_p=NMRapod;
[par2d,name2d,varseqFL]=tell2d(0);
zoomFL=0;

zinc=[8,1,-1,-8]; % zone increments
%button names 
guiname={'crs#1','crs#2','crs#3','Zoom','Unzoom','+8 zone','+1 zone','-1 zone','-8 zone',...
    'Apod.','Done'}; 

%initialize plot
cix=1;
nn=1:NMRseq(1,cix);
crs_n=[1,max(min(NMRseq(6,cix),nn(end)),1),nn(end)*.5]; %init curs. in n units
if (tsca),
  xx = nn*NMRseq(2,cix)*1e6; 
  crs=1e6*NMRseq(2,cix)*crs_n;
else
  xx = nn;
  crs=crs_n;
end
%apodization window (normalized)
apw=eval([apod_p.typ,'(NMRseq(1,cix),crs_n(2),1/(NMRseq(2,cix)*apod_p.fwid),'...
    'apod_p.half)']).*(crs_n(1)<=nn);
%cursor position matrix
xpos=crs'*[1 1];

%first plot		
funcs={'real','imag','abs'};
pid=-1*ones(1,3);
max_y=0;
clrs=['r','g','b'];
flg=1;
for ii=3:-1:1;
  if (ria_fl(ii)),
    yy=eval(sprintf('%s(TDfids(nn,cix))',funcs{ii}));
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
epid=zeros(1,5); %handles of apodwin & cursors: initialize
epid(4)=plot(xx,apw*max_y,'k');
epid(5)=line([xx(1) xx(end)],zeros(1,2),'color','k','linestyle','--');
if (tsca),
  xlabel('time (\mus)');
else
  xlabel('channels');
end
ylabel('amplitude (a.u.)');
ax=axis; ax(1:2)=[xx(1) xx(end)]; axis(ax);
epid(1)=line(xpos(1,:),ax(3:4),'color','g');
epid(2)=line(xpos(2,:),ax(3:4),'color','b');
epid(3)=line(xpos(3,:),ax(3:4),'color','m');
hold off
titid=title(sprintf('Zone %d, %s = %g',cix,name2d,par2d(cix)),...
  'fontsize',12);
set(gcf,'userdata',0); 

pos=[.08,.0,.08,.06]; 
dxpos=.14;
dypos=.115;
dpos=[dxpos,0];
% adding buttons to graph window, with dummy callbacks dummy1_0='dummy()' ...
for ii=1:prod(size(guiname));
  clbck=sprintf('set(gcf,''userdata'',%g);uiresume(gcf)',ii);  
  ui=uicontrol('style','push','units','normal','pos',pos,'string', ...
    guiname{ii},'call',clbck,'Tag','UINMR');
  if (ii==3),
    dpos=[0,-dypos];
    pos(1:2)=[.915,pos(2)+(length(guiname)-3)*dypos];
  end
  pos(1:2)= pos(1:2)+dpos;
end

figure(gcf);
%loop until Done button pressed (async polling)  
while(1),
  reapodFL=0; %1: rescale only, 2:recalc ydata; 3: xdata & ydata
  uiwait(gcf);  
  %callback definitions for gui buttons
  jj=get(gcf,'userdata');
  if (jj < 4), % moving a cursor 
    set(epid(jj),'linestyle','--');
    mpos=ginput(1);
    %validate position
    crs(jj)=max(xx(1),mpos(1));
    if (jj < 3) %3rd crs may be out of right bounds (no blc)
      crs(jj)=min(crs(jj),xx(end));
    end
    xpos(jj,:)=crs(jj)*[1 1];    
    set(epid(jj),'linestyle','-','xdata',xpos(jj,:));
    if (jj<3),
      reapodFL=2;
    end
  elseif (jj == 4), %zoom: use preblank and baseline cursors as wimdow limits
    %% zoomFl=1; 
    z=axis; z(1:2)=sort(crs([1 3])); axis(z);
  elseif (jj == 5), %unzoom
    %% zoomFL=0;
    z=axis; z(1:2)=[xx(1) xx(end)]; axis(z);
  elseif (jj < 10), %incrementing zone
    cix=max(min(cix+zinc(jj-5),size(TDfids,2)),1);
    if (varseqFL)
      reapodFL=3;
      nn=1:NMRseq(1,cix);
      if (tsca),
        xx=nn*NMRseq(2,cix)*1e6; 
      else
        xx=nn;
      end
      for k=1:3;
        crs(k)=min(max(xx(1),crs(k)),xx(end)); %validate cursors 
        xpos(k,:)=crs(k)*[1 1];                %(if out of bounds)
        set(epid(k),'xdata',xpos(k,:),'ydata',ax(3:4));
      end
      [ax,max_y]=replot(pid,epid,xx,cix);
      set(epid(5),'xdata',[xx(1) xx(end)]);
    else
      reapodFL=1;
      [ax,max_y]=replot(pid,epid,[],cix);
      for k=1:3;
        set(epid(k),'ydata',ax(3:4));
      end
    end
    set(epid(1:3),'visible','on');
    titid=title(sprintf('Zone %d, %s = %g',cix,name2d,par2d(cix)),...
      'fontsize',12);
  elseif (jj==10)
    reapodFL=2;
    apod_p=apodselect(apod_p);
  elseif (jj==11),
    % Done
    break
  end
  %recalculate and replot apodization window
  if (reapodFL),
    if (tsca), 
      crs_n=crs/(NMRseq(2,cix)*1e6);
    else
      crs_n=crs;
    end
    if (reapodFL >1),
      apw=eval([apod_p.typ,'(NMRseq(1,cix),crs_n(2),',...
          '1/(NMRseq(2,cix)*apod_p.fwid),apod_p.half)']).*(crs_n(1)<=nn);
    end
    if (reapodFL>2),
      set(epid(4),'xdata',xx,'ydata',apw*max_y,'visible','on');
    else
      set(epid(4),'ydata',apw*max_y,'visible','on');
    end
    
  end
 end

%clean-up: delete buttons
delete(findobj('Tag','UINMR'));

% restore scale in seconds instead of us.
if (tsca)
  crs=crs*1e-6;
end

%local macros:


%(private) REPLOT: replots fid from a new zone

function [ax,max_y]=replot(pid,epid,xvec,zone);
%pid=[id_real,id_imag,id_abs];
%epid=[handles of cursors, etc.];
%if xvec==[], xdata kept 
%zone: new zone index (no check whether it exists)
global TDfids NMRseq;

set(epid(1:4),'visible','off'); %I don't want these to affect scaling
N=NMRseq(1,zone);
fcn={'real','imag','abs'};
max_y=0;
fl=1;
%seek max_abs
for ii=3:-1:1;
  if (pid(ii)>0)
    y=eval(sprintf('%s(TDfids(1:N,zone))',fcn{ii}));
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




