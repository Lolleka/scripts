function crs=xfftsel(ria_fl,fsca,sh_fcn_ix); 
%% arguments:
%%   ria_fl=[real_flag,imag_flag,abs_flag] 
%%   fsca=1 -> uses time scale in Hz instead of channels
%%   sh_fcn_ix=1: shade real, =2: shade imag; =3: shade abs
%% output:
%%   crs=[center, width]

global FDfids NMRseq;

%% validation of arguments, apply defaults
if (nargin<1),
	ria_fl=[1 1 1];
elseif (size(ria_fl(:),1) < 3),
	ria_fl=[1 1 1];
elseif (~any(ria_fl))
	ria_fl=[1 1 1];
end;
if (nargin < 2),
  fsca=0;
end
autoshad=0; 
if (nargin > 2),
  sh_fcn_ix=max(1,min(sh_fcn_ix(1),3));
  if (~ria_fl(sh_fcn_ix)), %%inconsistent shading: switch to autodetect
    autoshad=1;
  end
else
  autoshad=1;
end
if (autoshad), %% auto-detect which curve is to be shaded
  xx=[2 1 3];
  [ii sh_fcn_ix]=max(xx.*(ria_fl~=0)); %% priority: abs>real>imag
end


%% program setup
[par2d,name2d,varseqFL]=tell2d(1);
zinc=[8,1,-1,-8]; %% zone incremnts
zoomx=0; 
%%button names 
guiname={'crs#1','crs#2','Zoom+','Zoom-','Unzoom','reset','+8 zone','+1 zone','-1 zone',...
	'-8 zone','Done'}; 

%%initialize plot
cix=1;
six=NMRseq(7,cix)/2 + 1;
nn=1:NMRseq(7,cix);
if (fsca),  
  xx = (nn - nn(end)/2 - 1)*NMRseq(8,cix); 
  d_x=NMRseq(8,cix);
  x0=d_x*(nn(end)/2+1);
else
  xx = nn - nn(end)/2 - 1;
  x0 = nn(end)/2 + 1;
  d_x=1;
end
initcrs=[0 0];
crs=initcrs;
xpos=[crs(1);crs(2)+crs(1)]*[1 1];


%first plot		
funcs={'real','imag','abs'};
pid=-1*ones(1,3);
epid=zeros(1,4); %handles of shading & cursors: initialize
clrs=['r','g','b'];
for ii=3:-1:1;
  if (ria_fl(ii)),
    yy=eval(sprintf('%s(FDfids(nn,cix))',funcs{ii}));
    pid(ii)=plot(xx,yy,clrs(ii));
    hold on;
    if (ii==sh_fcn_ix),
      xs=xx([six(1) six six(end)]);
      ys=[0 yy(six) 0];
      epid(3)=patch(xs,ys,clrs(sh_fcn_ix),'linestyle','none');
    end
  end
end
if (fsca),
  xlabel('Frequency (Hz)');
else
  xlabel('channels');
end
ylabel('amplitude (a.u.)');
ax=axis; ax(1:2)=[xx(1) xx(end)]; axis(ax);
epid(1)=line(xpos(1,:),ax(3:4),'color','b');
epid(2)=line(xpos(2,:),ax(3:4),'color','c');
epid(4)=line([xx(1) xx(end)],zeros(1,2),'color','k','linestyle','--');
%%%%%epid(3)=line(xpos(3,:),ax(3:4),'color','m');
hold off
titid=title(sprintf('Zone %d, %s = %g',cix,name2d,par2d(cix)),...
  'fontsize',12);
set(gcf,'userdata',0); 

%% install GUIs
pos=[.08,.0,.08,.06]; 
dxpos=.14;
dypos=.115;
dpos=[dxpos,0];
% add buttons to graphic window
for ii=1:prod(size(guiname));
  clbck=sprintf('set(gcf,''userdata'',%g);uiresume(gcf)',ii);  
  ui=uicontrol('style','push','units','normal','pos',pos,'string', ...
    guiname{ii},'call',clbck,'Tag','UINMR');
  if (ii==2),
    dpos=[0,-dypos];
    pos(1:2)=[.915,pos(2)+(length(guiname)-2)*dypos];
  end
  pos(1:2)= pos(1:2)+dpos;
end

figure(gcf);
%loop until Done button pressed (async polling)  
while(1),
  uiwait(gcf);  
  %callback definitions for gui buttons
  jj=get(gcf,'userdata');
  if (jj < 3), %% moving a cursor 
    set(epid(jj),'linestyle','--');
    mpos=ginput(1);
    xpos(jj,:)=min(max(xx(1),mpos(1)),xx(end))*[1 1];
    set(epid(jj),'linestyle','-','xdata',xpos(jj,:));
    if (jj==1),
      crs(1)=xpos(1,1);
      xpos(2,:)=min(max(xx(1),sum(crs)),xx(end))*[1 1];
      %%%ascis(sum(crs),cix); 
      set(epid(2),'xdata',xpos(2,:));
    end
    crs(2)=xpos(2,1)-crs(1);
  elseif (jj < 5), %zoom+-
    zoomx=zoomx-sign(jj-3.5);
    zoomx=max(0,min(round(log(NMRseq(7,cix))/log(2)-2),zoomx));
    zooma([xx(1) xx(end)],crs,zoomx);
  elseif (jj < 7), 
    zoomx=0; %unzoom
    axbox=axis;
    axis([xx(1) xx(end) axbox(3:4)])
    if (jj == 6),
      crs=initcrs;
      xpos=[crs(1);crs(2)+crs(1)]*[1 1];
      set(epid(1),'xdata',xpos(1,:));
      set(epid(2),'xdata',xpos(2,:));
    end
  elseif (jj < 11), %%incrementing zone
    cix=max(min(cix+zinc(jj-6),size(FDfids,2)),1);
    nn=1:NMRseq(7,cix);
    if (varseqFL)  %%sequence varies throughout zones: re-initializing
      if (fsca),  
        xx = (nn - nn(end)/2 - 1)*NMRseq(8,cix); 
        %%initcrs=[0 0];
        d_x=NMRseq(8,cix);
        x0=d_x*(nn(end)/2+1);
      else
        xx = nn - nn(end)/2 - 1;
        x0 = nn(end)/2 + 1;
        %%initcrs=[nn(end)/2+1 ,0];
        d_x=1;
      end
      %% validate cursors (if out of bounds)
      crs(1)=min(max(xx(1),crs(1)),xx(end)); 
      xpos(1,:)=crs(1)*[1 1];
      xpos(2,:)=min(max(xx(1),sum(crs)),xx(end))*[1 1];
      crs(2)=xpos(2,1)-crs(1);
      ax=replotft(pid,epid,xx,cix);
      for k=1:2;
        set(epid(k),'xdata',xpos(k,:),'ydata',ax(3:4));
      end
      set(epid(4),'xdata',[xx(1) xx(end)]);
  else
      ax=replotft(pid,epid,[],cix);
      for k=1:2;
        set(epid(k),'ydata',ax(3:4));
      end   
    end
    set(epid(1:2),'visible','on');
    titid=title(sprintf('Zone %d, %s = %g',cix,name2d,par2d(cix)),...
      'fontsize',12);
    if (zoomx)
      zooma([xx(1) xx(end)],crs,zoomx);
    end
  elseif (jj==11),
    %% Done
    break
  end
  
  if (jj<3 | jj==6 | (varseqFL)),  %%shading
    six=max(1,ceil((crs(1)-abs(crs(2))+x0)/d_x)):...
      min(nn(end),floor((crs(1)+abs(crs(2))+x0)/d_x));
    if isempty(six),
      six=max(1,min(nn(end),round((crs(1)+x0)/d_x)));
    end
  end
  reshade(epid(3),pid(sh_fcn_ix),six);
end

%clean-up: delete buttons
delete(findobj('Tag','UINMR'));
%validate result
crs(2)=abs(crs(2));
%%


%local macros:

%(private) REPLOTFT: replots fid from a new zone
function ax=replotft(pid,epid,xvec,zone);
%pid=[id_real,id_imag,id_abs];
%epid=[handles of cursors, etc.];
%if xvec==[], xdata kept 
%zone: new zone index (no check whether it exists)
global FDfids NMRseq;

set(epid(1:3),'visible','off'); %I don't want these to affect scaling
N=NMRseq(7,zone);
fcn={'real','imag','abs'};
fl=1;
for ii=3:-1:1;
  if (pid(ii)>0)
    y=eval(sprintf('%s(FDfids(1:N,zone))',fcn{ii}));
    if (length(xvec)),
      set(pid(ii),'Xdata',xvec,'Ydata',y);
    else
      set(pid(ii),'Ydata',y);
      if (fl),
        xvec=get(pid(ii),'Xdata');
        fl=0;
      end
    end
  end
end
axis auto;
ax=axis; ax(1:2)=[xvec(1),xvec(end)]; axis(ax);


%(private) : RESHADE:  redraw shading 
function reshade(shid,plid,ix);
x=get(plid,'xdata');
y=get(plid,'ydata');
xs=x([ix(1) ix ix(end)]);
ys=[0 y(ix) 0];
set(shid,'xdata',xs,'ydata',ys,'visible','on');


%(private) : ZOOMA
function zooma(x_ax,crs,zoomx)
axbox=axis;
axbox(1:2)=[min(max(x_ax(1),crs(1)-(crs(1)-x_ax(1))/2^zoomx),crs(1)-crs(2)) ...
    max(min(x_ax(2),crs(1)+(x_ax(2)-crs(1))/2^zoomx),crs(1)+crs(2))];
axis(axbox);
