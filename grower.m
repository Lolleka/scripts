%    N P K
teasp = 4.92892159; %ml
ref_two_lits = 2/3.79;

FM= [5 0 1];
FG= [2 1 6];
FB= [0 5 4];

CV = [7.5 3 6];
CF = [3 4 5];

GHplan = [2 1 1; 5 4 1; 5 4.5 1.5; 5 5 2; 4 4 4; 4 1 5; 4 1 5; 4 1 6; 4 1 6; 4 1 6; 4 0 8; 4 0 8; 2.5 0 8]
Weeks = 1:size(GHplan,1);
GHStrength = sum(GHplan,2);
Nitrogen = GHplan*[FM(1) FG(1) FB(1)]' ./ GHStrength;
Phosphorus = GHplan*[FM(2) FG(2) FB(2)]' ./ GHStrength;
Potassium = GHplan*[FM(3) FG(3) FB(3)]' ./ GHStrength;

 
figure(1);
plot(Weeks,Nitrogen,'red');
hold on;
plot(Weeks,Phosphorus,'green');
plot(Weeks,Potassium,'blue');
plot(Weeks, GHStrength,'black');
hold off;

%in this case the third number is plain water+
%Cplan = [ 1.8 0.2 2; 2 0 2];
Cplan=[];
[VV,VF,VW]= meshgrid(0:60,0:60,0:60);
[N, P, K] = Compo(VV,VF,VW);
for i = 1:size(Weeks,2)
TN = Nitrogen(i);
TP = Phosphorus(i);
TK = Potassium(i);
npkdist = sqrt((N - TN).^2 + (P-TP).^2 + (K - TK).^2);
[opt, optI] = min(npkdist(:));
[Ia, Ib, Ic] = ind2sub(size(npkdist),optI)

Cplan = [Cplan; [VV(Ia,Ib,Ic), VF(Ia,Ib,Ic), VW(Ia,Ib,Ic)]];

end

Weeks = 1:size(Cplan,1);
CStrength = sum(Cplan,2);
Nitrogen = Cplan*[CV(1) CF(1) 0]' ./ CStrength;
Phosphorus = Cplan*[CF(2) CF(2) 0]' ./ CStrength;
Potassium = Cplan*[CV(3) CF(3) 0]' ./ CStrength;
figure(1);
hold on;
plot(Weeks,Nitrogen,'.','color','red');
plot(Weeks,Phosphorus,'.','color','green');
plot(Weeks,Potassium,'.','color','blue');
%plot(Weeks, CStrength,'.','color','black');
hold off;

Cplan = [Cplan(:,1).*GHStrength./CStrength ...
         Cplan(:,2).*GHStrength./CStrength ...
         Cplan(:,3).*GHStrength./CStrength]

Cplan_2l_teasp = chop(Cplan *ref_two_lits/ teasp,1)
chop(Cplan_2l_teasp,1)*teasp
%figure(2);
