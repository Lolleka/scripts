function config;
instdir=which('setup');
if (length(instdir)),
  N=max(find(instdir=='/'));
  if (~isunix),
    N=max([N find(instdir=='\')]);
  end
  instdir=instdir(1:N);
else
  error('The NMR2001 package is not in your Matlab path');
end

vnames={'NMRpath' 'NMRplotsty' 'NMRapod' 'NMRfftwin' 'NMRfitflags' 'NMRstat'};

for k=1:length(vnames);
  eval(['global ',vnames{k}]);
end

if (~isstruct(NMRplotsty)),
  error('Uninitialized workspace. Run "setup" before');
end

markers={'.' 'o' 'x' '+' '*' 's' 'd' 'v' '^' '<' '>' 'p' 'h'};
lines={'-' ':' '-.' '--';'solid' 'dotted' 'dash-dotted' 'dashed'};
colors={' ' 'y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'; 'default' 'yellow' 'magenta' ...
	'cyan' 'red' 'green' 'blue' 'white' 'black'};


udt=struct('id',zeros(1,16),'sel',zeros(1,16),'vl',zeros(1,16),...
	   'path',NMRpath);
txth=zeros(1,11); %text-type objects
fig_id=figure;
set(fig_id,'position',[90,90,600,450])
if isunix,
  fs=11;
else
  fs=8;
end

%%% CALLBACKS
%checkbuttons, popup menus
clbck=['N=get(gcbo,''userdata'');val=get(gcbo,''value'');u=get(gcf,' ...
	'''userdata'');u.sel(N)=1;u.vl(N)=val;set(gcf,''userdata'',u)'];

%apodization selection
clbck2=['[val st]=apodselect([],1);N=get(gcbo,''userdata'');u=' ...
	' get(gcf,''userdata'');u.sel(N)=st;set(gcf,''userdata'',u)'];

%editable numerics
clbck3 = ['N=get(gcbo,''userdata'');u=get(gcf,''userdata'');',...
	  'r=sscanf(get(gcbo,''string''),''%f'');',...
	  'if(isstr(r)), set(gcbo,''string'',num2str(u.vl(N)));',...
	  'else, u.vl(N)=r(1);u.sel(N)=1;set(gcf,''userdata'',u);end'];

%editable path
clbck4= ['u=get(gcf,''userdata'');r=get(gcbo,''string'');if (length(r)),',...
	 'val=(exist(r)~=7);else val=0;end; if (val), set(gcbo,',...
	 '''string'',u.path);else, N=get(gcbo,''userdata'');u.path=r;',...
	 'u.sel(N)=1;set(gcf,''userdata'',u);end'];


x0=.04; dx=.125;
y0=.78; dy=.125;

%%% TEXT LABELS
%% Main label
N=1;
txth(1)=uicontrol('style','text','units','norm','string',...
		  'Configuration Panel', 'Position',[.3,.93,.4,.05], ...
		  'fontsize', fs,'Foreground','b','Horiz','Center');

%% Red labels
labs={'Plot Setup','Fitting setup','Apodization',['Default FFT window'],...
      'Data Path'};
pos=ones(length(labs),1)*[.05,.88,.15,.05]; 
pos(2:end,2)=pos(2:end,2)-dy*1.1*(2.415+(0:size(pos,1)-2))';
pos(4,4)=pos(4,4)*1.5;
for ii=1:length(labs);
  N=N+1;  
  txth(N)=uicontrol('style','text','units','norm','string',...
		    labs{ii}, 'Position',pos(ii,:), ...
		    'fontsize', fs,'Foreground','r','Horiz','Center');
end

%% Black labels
pos=[x0              y0 .17 dy*2/3
     x0+.185+2.15*dx y0 .13 dy*2/3
     ones(3,1)*[x0+1.4*dx y0-dy*3.95 dx*.56 dy/3]];
pos(4:end,1)=pos(4:end,1)+dx*1.9*(1:2)';
labs={{'Marker symbols';'& colors'},{'Line style';'& color'}, ...
      'Shift','Width','Units'};
for ii=1:length(labs);
  N=N+1;  
  txth(N)=uicontrol('style','text','units','norm','string',labs{ii},...
		    'Position',pos(ii,:),'fontsize',fs,'Horiz','Center');
end

%%% POPUP GUIs
%% Plot Setup: markers & line style

M=0;
x0=x0+.185; 
pos=ones(4,1)*[x0 y0+dy/16 dx*0.525 dy/2];
pos(2:4,[1 3])=[x0+dx*.65 dx*1.15; x0+dx*2.15+.145 dx*1.35; ...
		x0+dx*3.625+.145 dx*1.15];
labs={'markers','colors','lines','colors'};
substring={strtok(NMRplotsty.sym,cat(2,colors{1,:}))
	   strtok(NMRplotsty.sym,cat(2,markers{:}))
	   strtok(NMRplotsty.line,cat(2,colors{1,:}))
	   strtok(NMRplotsty.line,cat(2,lines{1,:})) };

for ii=1:length(labs);
  M=M+1;
  val=[find(strcmp(substring{ii},eval([labs{ii},'(1,:)']))) 1];
  udt.vl(M)=val(1);
  udt.id(M)=uicontrol('style','popup','units','norm','position',...
		      pos(ii,:),'value',val(1),'string',...
		      eval([labs{ii},'(end,:)']),...
		      'fontsize',fs+1,'userdata',M,'callback',clbck);
end

%%% CHECK BUTTONS
%% Plot & Fitting Flags
x0=.04; 
y0=y0-dy*.85;
pos=ones(6,1)*[x0 y0+dy/16 dx*1.15 dy/2];
pos([2 4:6],3)=dx*[1.8;1.75;2.25;2.5];
pos(2:end,1)=x0+dx*[1.5;3.65;5.15;1.58;4.04];
pos(5:end,2)=pos(1,2)-dy*1.125*[1; 1];

val=[(NMRplotsty.ebar~=0) (NMRplotsty.drawpeaks~=0) (NMRplotsty.grid~=0) ...
     (NMRplotsty.ylog~=0) (NMRfitflags.batch~=0) (NMRfitflags.offs~=0)];
labs={'Error bars', 'Draw Components', 'Draw Grid', 'T1/T2: Y Log Scale',...
      'Automatic fit (batch mode)', 'Enable offset components'};

for ii=1:length(labs);
  M=M+1;
  udt.vl(M)=val(ii);
  udt.id(M)=uicontrol('style','checkbox','units','norm','position',...
		      pos(ii,:),'value',val(ii),'string',labs{ii},...
		      'fontsize',fs,'userdata',M,'callback',clbck);
end

%%% PUSH BUTTONS
%% Apodization
y0=y0-2.05*dy;
x0=.05+1.75*dx; 
M=M+1;		
udt.id(M)=uicontrol('style','push','units','norm','position',...
		    [x0 y0-dy/16 dx*3 dy/2],'string',...
		    'Change apodization window','fontsize',fs, ...
		     'userdata',M,'callback',clbck2);

%%% EDITABLE TEXT
%%FFT window (numeric strings)
y0=y0-1.09*dy;
x0=.05+1.92*dx; 
pos=ones(2,1)*[x0 y0 dx*1.125 dy*.4];
pos(2,1)=pos(2,1)+dx*1.9;
val=[NMRfftwin.foffs(1), NMRfftwin.fwid(1)];
for ii=1:length(val);
  M=M+1;
  udt.vl(M)=val(ii); 
  udt.id(M)=uicontrol('style','edit','units','norm','string',...
		      num2str(val(ii)),'call',clbck3,'position',pos(ii,:), ...
		      'Horiz','center','fontsize',fs,'userdata',M);
end
%path
pos=[x0-dx/2 y0-1.182*dy dx*4 dy*.4];
M=M+1;		
udt.id(M)=uicontrol('style','edit','units','norm','string',...
		    NMRpath,'call',clbck4,'position',pos, ...
		    'Horiz','left','fontsize',fs,'userdata',M);

%specialized callback for interactive path browsing
clbck8 = ['u=get(gcf,''userdata'');[r,val]=uigetfile([u.path,''\*.*''],',...
	  '''NMRpath''); if(isstr(val)), u.sel(',int2str(M),')=1;',...
	  'u.path=val;set(u.id(',int2str(M),'),''string'',val);', ...
	  'set(gcf,''userdata'',u); end;'];


%BUTTONGROUP: frequency/channels units 
icon_opt=sprintf('''Horiz'',''center'',''fontsize'',%d',fs);

icona = ['text(.5,.5,''   Hz   '',',icon_opt,')'
	 'text(.5,.5,''Channels'',',icon_opt,')'
	];
inistat=zeros(1,2); 
inistat(2-(NMRfftwin.funit==0))=1;
pos=[x0+3.8*dx y0 dx*1.64 dy*.45];

M=M+1;
%buttongroup: specialized callback
clbck9 = ['u=get(gcf,''userdata'');tg=get(gcbo,''tag'');',...
	  'u.vl(',int2str(M),')=eval(tg(2)); u.sel(',int2str(M),')=1; ',...
		    'set(gcf,''userdata'',u)'];
udt.vl(M)=2-(NMRfftwin.funit==0); 
udt.id(M)=btngroup('GroupID', 'Funits', 'ButtonID', ...
		   ['U1';'U2'],'Callbacks',clbck9, ...
		   'exclusive','yes',... %%'Userdata',M,...
		   'groupsize',[1,2],'Position', pos, ...
		   'IconFunctions', icona,'InitialState',inistat);

		   
%%% PUSHBUTTON: browse path
x0=x0+dx*3.7;
y0=y0-1.184*dy;

M=M+1;
udt.id(M)=uicontrol('style','push','units','norm','position',...
		    [x0 y0 dx dy*.5],'string','Browse','fontsize',fs, ...
		     'userdata',M,'callback',clbck8);


%%% THE "DONE" BUTTON
done_id=uicontrol('style','push','units','norm','call','delete(gcbo)', ...
		  'string','Done','Position',[.1 .04 dx*1.15 dy*.55],...
		  'fontsize',fs+1);

set(fig_id,'userdata',udt);
waitfor(done_id);
udt=get(fig_id,'userdata');
if (~length(udt)),
  warning('Abnormal GUI exit');
  close(fig_id);
  return
end

NMRplotsty.sym=[markers{udt.vl(1)} colors{1,udt.vl(2)}];
NMRplotsty.line=[lines{1,udt.vl(3)} colors{1,udt.vl(4)}];
NMRplotsty.ebar=udt.vl(5);
NMRplotsty.drawpeaks=udt.vl(6);
NMRplotsty.grid=udt.vl(7);
NMRplotsty.ylog=udt.vl(8);

NMRfitflags.batch=udt.vl(9);
NMRfitflags.offs=udt.vl(10);

NMRfftwin.foffs=udt.vl(12);
NMRfftwin.fwid=udt.vl(13);
NMRfftwin.funit=udt.vl(15)-1;

NMRpath=udt.path;
close(fig_id);
if (any(udt.sel)),
  val=questdlg('Save configuration?','NMR2001 Config','Yes','No','Yes');
  if (strcmp(val,'Yes')),
    if (isunix),
      [M,N]=unix(['touch ',instdir,'config.mat']);
      if (M),
	disp(['You don''t have write permission in ',instdir,...
		 ', trying with current directory']);
	instdir=pwd;
      end
    end
    save([instdir,'/config.mat'],vnames{:});
    disp(['Setup saved into ',instdir,'/config.mat']);
  end  
end

%% Copyright (C) G. Allodi
%% Revision 2.0   14-sept-2002
