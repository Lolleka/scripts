function [ria_,sca_,fd_]=plotmodesel_(ria,sca,fd)
sca=(sca~=0);
fd=(fd~=0);
ria=(ria~=0);
if (all(ria==0)),
  ria(1+2*fd)=1;
end
ria_=ria;
sca_=sca;
fd_=fd;


udt=struct('flg',[0 (1+(sca~=0)) (1+(fd~=0))],'ria',ria,'id',zeros(1,6));

set(figure,'position',[100,100,400,300]);
if isunix,
  fs=12;
else
  fs=9;
end

tid=uicontrol('style','text','units','norm','string','Display setup',...
	      'Position',[.3,.93, .4,.06],'fontsize', fs,'Foreground','b',...
	      'Horiz','Center');


x0=.3; dx=.18;
y0=.74; dy=.15;

%% callback definition

% calback for radio-button-style groups (exclusive)
clbck = ['u=get(gcf,''userdata'');tg=get(gcbo,''tag'');',...
	 'u.flg(abs(tg(1))-64)=eval(tg(2));set(gcf,''userdata'',u)'];


%callback for Real, ... Abs, selection (non-exclusive)
clbck1 = ['u=get(gcf,''userdata''); u.ria=btnstate(gcf,''Tracks'');',...
	  'if all(u.ria==0), df=1+2*(u.flg(2)>1); u.ria(df)=1;',...
	  'btnstate(''set'',gcf,''Tracks'',u.ria); end;',...
	   'set(gcf,''userdata'',u)'];

icon_opt=sprintf('''Horiz'',''center'',''fontsize'',%d',fs);
icon1 = [
    'text(.5,.5,''Real'',',icon_opt,')'
    'text(.5,.5,''Imag'',',icon_opt,')'
    'text(.5,.5,''Abs '',',icon_opt,')'
        ];

icon2 = [
    'text(.5,.5,''channels'',',icon_opt,')'
    'text(.5,.5,''us / MHz'',',icon_opt,')'
        ];

icon3 = [
    'text(.5,.5,''  Time   '',',icon_opt,')'
    'text(.5,.5,''Frequency'',',icon_opt,')'
        ];


udt.id(1)=btngroup('GroupID', 'Tracks', 'ButtonID', ...
                 ['A1';'A2';'A3'],'Callbacks',clbck1, ...
                 'groupsize',[1,3],'Position', [x0 y0 dx*3 dy ],...
                 'IconFunctions', icon1,'InitialState',ria ~= 0);

udt.id(4)=uicontrol('style','text','units','norm','string','Tracks',...
		    'Position',[.04 y0+dy/4 .2 dy/2],'fontsize',fs,...
		    'Horiz','Center');

y0=y0-dy*1.33;

udt.id(2)=btngroup('GroupID', 'XUnits', 'ButtonID', ...
                 ['B1';'B2'],'Callbacks',clbck, ...
		 'exclusive','yes',...
                 'groupsize',[1,2],'Position', [x0 y0 dx*2.5 dy ],...
                 'IconFunctions', icon2,'InitialState',rem([1 0] + sca,2));
       
udt.id(5)=uicontrol('style','text','units','norm','string','X-axis units',...
		    'Position',[.04 y0+dy/4 .2 dy/2],'fontsize',...
		    fs,'Horiz','Center');

y0=y0-dy*1.33;

udt.id(3)=btngroup('GroupID', 'Domain', 'ButtonID', ...
                 ['C1';'C2'],'Callbacks',clbck, ...
		 'exclusive','yes',...
                 'groupsize',[1,2],'Position', [x0 y0 dx*2.5 dy ],...
                 'IconFunctions', icon3,'InitialState',rem([1 0] + fd,2));

udt.id(6)=uicontrol('style','text','units','norm','string',{'Time/Freq. ',...
		    'domain'},'Position',[.04 y0+dy/6 .2 dy*.7],'fontsize',...
		    fs,'Horiz','Center');


y0=y0-dy*1.33;	  

done_id=uicontrol('style','push','units','norm','call','delete(gcbo)', ...
		  'string','Done','Position',[x0+2*dx y0 dx dy]);

set(gcf,'userdata',udt);
waitfor(done_id);
usrdt=get(gcf,'userdata');       

if (~length(usrdt)),
  warning('Abnormal GUI exit, applying defaults');
  close(gcf);
  return
end
ria_=usrdt.ria;
sca_=(usrdt.flg(2) > 1);
fd_=(usrdt.flg(3) > 1);
close(gcf);

%private function, changes plot options interactively from within 
%nmrdisplay.  