function [ap,status]=apodselect(ap_i,globcfg);
global NMRapod %NMRapod = struct('typ','hygausswin','fwid',1e4,'half',0)

if (~exist('ap_i')),
  ap_i=NMRapod;
end
if (~isstruct(ap_i)), %not a struct 
  ap_i=NMRapod;
end;
if (nargin < 2), 
  globcfg=0; %flag: if == 1, global setup, disable query 
end
if (length(fieldnames(ap_i)) ~= length(fieldnames(NMRapod))), %invalid struct
  ap_i=NMRapod;
end;

%known window types 
types={'expowin','gausswin','hygausswin','rectwin','flatwin'};
tyx=find(strcmp(ap_i.typ,types));
if (~length(tyx)),
  tyx=5;
end
tyx=tyx(1);
fhx=1+ap_i.half;


%%Matlab specific code
udt=struct('id',[],'flg',[tyx ((ap_i.half~=0)+1) 1 ap_i.fwid],'st',0);
set(figure,'position',[100,100,480,360])
if isunix,
  fs=11;
else
  fs=8;
end

udt.id=[udt.id,uicontrol('style','text','units','norm','string',...
			 'Apodization setup', 'Position',[.3,.93, ...
		    .4,.06],'fontsize', fs,'Foreground','b','Horiz','Center')];

x0=.25; dx=.15;
y0=.76; dy=.125;
icon_opt=sprintf('''Horiz'',''center'',''fontsize'',%d',fs);

icon1 = [
    'text(.5,.5,''  Expo   '',',icon_opt,')'
    'text(.5,.5,''  Gauss  '',',icon_opt,')'
    'text(.5,.5,'' HypGauss'',',icon_opt,')'
    'text(.5,.5,''Rectangle'',',icon_opt,')'
    'text(.5,.5,''  None   '',',icon_opt,')'
	];

icon2 = ['text(.5,.5,''Full'',',icon_opt,')'
	 'text(.5,.5,''Half'',',icon_opt,')'
	];

icon3 = ['text(.5,.5,''Local '',',icon_opt,')'
	 'text(.5,.5,''Global'',',icon_opt,')'
	];

inbut1=zeros(1,5); inbut1(tyx)=1;
inbut2=zeros(1,2); inbut2(udt.flg(2))=1;
inbut3=zeros(1,2); inbut3(1)=1;
clbck = ['u=get(gcf,''userdata'');tg=get(gcbo,''tag'');',...
	 'u.flg(abs(tg(1))-64)=eval(tg(2));u.st=1;set(gcf,''userdata'',u)'];

clbck4 = ['u=get(gcf,''userdata'');r=sscanf(get(gcbo,''string''),''%f'');',...
	  'if(isstr(r)), set(gcbo,''string'',num2str(u.flg(4)));',...
	  'else, u.flg(4)=r(1);u.st=1;set(gcf,''userdata'',u);end'];

udt.id=[udt.id,...
	btngroup('GroupID', 'WinType', 'ButtonID', ...
		 ['A1';'A2';'A3';'A4';'A5'],'Callbacks',clbck, ...
		 'exclusive','yes',...
		 'groupsize',[1,5],'Position', [x0 y0 dx*5 dy ],...
		 'IconFunctions', icon1,'InitialState',inbut1)
       ];

udt.id=[udt.id,uicontrol('style','text','units','norm','string',['Window' ...
		    ' type'],'Position',[.04 y0+dy/4 .2 dy/2],'fontsize', ...
			 fs,'Horiz','Center')
       ];

y0=y0-dy*1.33;

udt.id=[udt.id,...
	btngroup('GroupID', 'IsHalf', 'ButtonID', ...
		 ['B1';'B2'],'Callbacks',clbck, ...
		 'exclusive','yes',...
		 'groupsize',[1,2],'Position', [x0 y0 dx*2 dy ],...
		 'IconFunctions', icon2,'InitialState',inbut2)
       ];

udt.id=[udt.id,uicontrol('style','text','units','norm','string',['Full/' ...
		    'Half window'],'Position',[.04 y0+dy/5 .2 dy*.8],...
			 'fontsize',fs,'Horiz','Center')
       ];

y0=y0-dy*1.33;
if (~globcfg), %%disabled with global setup
  udt.id=[udt.id,...
	  btngroup('GroupID', 'IsLocal', 'ButtonID', ...
		   ['C1';'C2'],'Callbacks',clbck, ...
		   'exclusive','yes',...
		   'groupsize',[1,2],'Position', [x0 y0 dx*2 dy ],...
		   'IconFunctions', icon3,'InitialState',inbut3)
	 ];

  udt.id=[udt.id,uicontrol('style','text','units','norm','string',...
			   'Setup Scope','Position',[.04 y0+dy/4 .2 dy/2],...
			   'fontsize',fs,'Horiz','Center')
	 ];
end
y0=y0-dy*1.33;
udt.id=[udt.id,...
	uicontrol('style','edit','units','norm','string',...
		  num2str(udt.flg(4)),'call',clbck4,'position',[x0+2*dx ...
		    y0 dx dy],'Horiz','center','fontsize',fs)];

udt.id=[udt.id,...
	uicontrol('style','text','units','norm','string',...
		  {'Window width [Hz]','(meaning depends on win type)'},...
		  'Position',[.04 y0+dy/5 .2+2*dx dy*.8],...
		  'fontsize',fs,'Horiz','Center')
       ];

y0=y0-dy*1.33;	  

done_id=uicontrol('style','push','units','norm','call','delete(gcbo)', ...
		  'string','Done','Position',[x0+2*dx y0 dx dy]);
set(gcf,'userdata',udt);

waitfor(done_id);

usrdt=get(gcf,'userdata');
if (~length(usrdt)),
  warning('Abnormal GUI exit, applying defaults');
  ap=ap_i;
  close(gcf);
  return
end
flg=usrdt.flg; %return value
if (nargout > 1),
  status = usrdt.st; %%status=0 if nothing has been changed 
end
close(gcf);

ap.typ=types{flg(1)};
ap.half=(flg(2)>1);
ap.fwid=flg(4);
if (globcfg | flg(3)>1),
  NMRapod=ap;
end


