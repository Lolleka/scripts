
function [res_x res] = ZFCFit(filename, TCut1, TCut2, nu0, nuT, renorm, x0)

%data = load('Prova.dat');
data = load(filename);
%TCut = data(find(data(:,2) == max(data(:,2))),1)*2/5;
%TCut = TCut(1);
%TCut = 70;
ind = find(data(:,1)> TCut1 & data(:,1)< TCut2);
fit_data = [data(ind,1) data(ind,2)];

F = @(x,xdata)(renorm)*CalcZFC(xdata,x(1),nu0,nuT,x(2),x(3));

%figure(1);
%plot (data(:,1),data(:,2));
%hold on;

[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,fit_data(:,1)',fit_data(:,2)');

%plot(data(:,1),F(x,data(:,1)'),'Color','red');

res_x = x;
res = [data(:,1)'; F(x,data(:,1)')];
%hold off;
%1E-30
%figure(1);
%for i=100:10:200
 %   i
 %   r = CalcZFC(1E-25,1E6,2,i,50);
 %   plot(r(1,:),r(2,:));
 %   hold on;
%end
%hold off;
end