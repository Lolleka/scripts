rho_IO(1) = 4900;  %Density of Iron Oxide (Kg/m^3)
rho_IO(2) = 4900;  %Density of Iron Oxide (Kg/m^3)
rho_O(1) = 1150;   %Density of organic component (Kg/m^3)
rho_O(2) = 1150;   %Density of organic component (Kg/m^3)
r(1) = 12.7E-9;
r(2) = 11.6E-9;
D(1) = 78.6E-9;
D(2) = 121.8E-9;
%Iron Fractions:

fm_IO(1) = 0.868;    % % of maghemite
fm_IO(2) = 0.921;    % % of maghemite

fv_IO = 1./((1-fm_IO).*rho_IO./(rho_O.*fm_IO)+1); %volume fraction of Iron Oxide
%N = round(fv_IO .* (D.^3) ./ (r.^3)); %Number of particles in a cluster

%Calculation of Mv

MS(1) = 71.2;   %emu/gr
MS(2) = 73.4;   %emu/gr

MPart = MS .* rho_IO;   % 1E-3 emu/cm^3
MV = MPart .* fv_IO; %.* fv_IO     % A/m

%Calculation of Delta_omega*tau_D
%hydrodynamic diameters
d_hyd(1) = 78.6E-9;  %m
d_hyd(2) = 121.8E-9; %m
%particle diameters:
d(1) = 50E-9;
d(2) = 86E-9;

gamma_P = 2 * pi * 42.576E6; %Hz*rad/T
mu_0 = 4*pi*1E-7; %T*m/A
DD = 3E-9; %m^2/s
Delta_omega = (gamma_P * mu_0) * MV ./ 3
tau_D = (d.^2) ./ (4*DD)

regime = Delta_omega.*tau_D
% Relaxivities
r2_exp(1) = 405.1;
r2_exp(2) = 508.3;

nu_mat = 1.57E-5;
r2_theo = 2*pi*gamma_P*mu_0*nu_mat*MV/(9*sqrt(3));

myres_exp = r2_exp.*fv_IO./(MV.^2);
myres_theo = r2_theo.*fv_IO./(MV.^2);

r2_theo./r2_exp;

%nu_mat = r2_exp*(9*sqrt(3))./(MV*2*pi*gamma_P*mu_0)
data_vuong = load('Data_Vuong.dat', '-ascii');
fit_vuong = load('Fit_Vuong.dat', '-ascii');
%Fe_conc = fv_IO./nu_mat

figure(1);
data_vuong_p = plot(data_vuong(:,1), data_vuong(:,2), 'o');
hold on;
fit_vuong_p = plot(fit_vuong(:,1), fit_vuong(:,2));
myres_exp_p = plot(d*1E9,myres_exp,'o');
%myres_theo_p = plot(d*1E9,myres_theo,'o');
hold off;
set(data_vuong_p,'MarkerEdgeColor','black','MarkerFaceColor','green','MarkerSize',3);
set(fit_vuong_p,'Color','black');
set(gca,'XScale','log','YScale','log');


ENDOREM_fv_IO = 0.23;
ENDOREM_MV = 0.77E5 .* ENDOREM_fv_IO;
ENDOREM_d = 80e-9;
ENDOREM_r2_exp = 107;
ENDOREM_tau_D = ENDOREM_d .^2 ./ DD;
ENDOREM_Delta_omega = (gamma_P * mu_0) * ENDOREM_MV ./ 3;

ENDOREM_regime = ENDOREM_Delta_omega.*ENDOREM_tau_D
ENDOREM_myres_exp = ENDOREM_r2_exp.*ENDOREM_fv_IO./(ENDOREM_MV.^2) %(0.77E5^2) %(ENDOREM_MV.^2)
hold on;
myres_exp_p = plot(ENDOREM_d*1E9,ENDOREM_myres_exp,'o','Color','red');
hold off;

other_fv_IO = [0.3 0.5];
other_MV = [1.23E5 1.79E5] .* other_fv_IO;
other_d = [34e-9 63e-9];
other_r2_exp = [540 630];
other_tau_D = other_d .^2 ./ DD;
other_Delta_omega = (gamma_P * mu_0) * other_MV ./ 3;

other_regime = other_Delta_omega.*other_tau_D
other_myres_exp = other_r2_exp.*other_fv_IO./(other_MV.^2)
%hold on;
%myres_exp_p = plot(other_d*1E9,other_myres_exp,'o','Color','red');
%hold off;