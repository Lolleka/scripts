%F = @(x) BPPDistribution(x, {300, 1E-9, 2, 0, 10, 0, 1E7})
T1vsT = load('sample.dat');
T = T1vsT(:,1);
T1 = T1vsT(:,2);
InvT1 = 1 ./ T1;

figure (1);
clf; % Clear the last figure (in this case figure 1)

plot(T,InvT1,'-o');


initial_tau0 = 2E-12;
initial_DPos = 130;
initial_DWidth = 5;
initial_Ampl = 7E8;
H = 0.5; %Tesla

t = [1E-12 1 1 1E8];
p0 = [initial_tau0 initial_DPos initial_DWidth initial_Ampl] ./ t;


omegaL = 42.576E6*2*pi*H;


F = @(p,xdata)BPPDistribution(xdata,p(1)*t(1),p(2)*t(2),p(3)*t(3),0,p(4)*t(4),omegaL);

T_subset = T(1:20);
InvT1_subset = InvT1(1:20);

[p,resnorm,~,exitflag,output] = lsqcurvefit(F,p0,T_subset,InvT1_subset);

p .* t

BestFit = BPPDistribution(T,p(1)*t(1),p(2)*t(2),p(3)*t(3),0,p(4)*t(4),omegaL);

hold on;
plot( T, BestFit, '-', 'Color', 'red');
