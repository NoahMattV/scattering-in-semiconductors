% Noah Van Der Weide

% GaAs Scattering Rates
% 1) Acoustic Phonon Scattering (abs and em together within elastic and 
% equipartition approx.)
% 
% 2) Polar optical phonon scattering (abs and em separately for electrons
% in the Gamma valley
%
% 3) Intervalley Scattering (G-L and L-G, together in one graph; 
% G-X and X-G together, X-L and L-X together, L-L, X-X)
%
% 4) Ionized Impurity Scattering at doping densities 10^17 and 10^19 cm^-3
%
% Part B) Plot Relaxation Rates and Momentum Relaxation Rates for POP and 
% Ionized Impurity Scattering in the same energy range 
% (to show how different they are)
%
% Electron energy range should be from 0 - 2 eV. 

close all;
clear;

e = 1.602e-19; % charge of an electron
q = e;
ko = 12.9; % low freq. dielectric const.
kinf = 10.92; % high freq. dielectric const.
%hbar = 6.582e-16; % Reduced Planck's Constant (eV*s)
hbar = 1.054e-34; % (J*s)
%rho = 5.36; % mass density (g/cm^3)
rho = 5360; % (kg/m^3)
a = 5.462e-10; % lattice constant (m)
kT = 0.0259*e; % (J)
%vs = 5.24e5; % Longitudinal acoustic velocity (cm/s)
vs = 5240; % (m/s)

Dac_G = 7.01*e; % Electron acoustic deformation potential (eV) - Gamma
Dac_L = 9.2*e;
Dac_X = 9.0*e;

mo = 9.11e-31; %kg
mdos = 0.067*mo;
Ec = 0;
E = linspace(0,2*e,1001); % 0 to 2 eV

%Z = 64; % Protons in GaAs
Z = 1;
Es = 12.9; % relative dielectric constant
Eo = 8.854e-12; %F/m vacuum permittivity
Ep = Es*Eo; 

% ---------------------------------------------------------------------------------------
% 1) Acoustic Phonon Scattering - within elastic and equipartition approximations
% ---------------------------------------------------------------------------------------

gc3d = zeros(1001);
GG = zeros(1001);
GL = zeros(1001);
GX = zeros(1001);
G = zeros(1001);

for i = 1:1001
   %gc3d(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * sqrt(E(i) - Ec); % 3D DOS for GaAs
   %G(i) = 2 * 2*pi/hbar * Dac_G^2 * kT/(2*rho*vs^2) * (1/2) * gc3d(i);
   GG(i) = (Dac_G^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*mdos/(hbar^2))^(3/2) * E(i)^(1/2);
   GL(i) = (Dac_L^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*mdos/(hbar^2))^(3/2) * E(i)^(1/2);
   GX(i) = (Dac_X^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*mdos/(hbar^2))^(3/2) * E(i)^(1/2);
   G(i) = GG(i) + GL(i) + GX(i);
end

figure();
plot(E,G(1:1001));
title('Acoustic Phonon Scattering - GaAs');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');

% ---------------------------------------------------------------------------------------
% 2) Polar Optical Phonon Scattering - Absorption and Emission separately
% for electrons in the Gamma valley
% ---------------------------------------------------------------------------------------
% Strongly inelastic (Ge != 0)
% Strongly anisotropic (Gm != G)

Gm_pop_abs = zeros(1001); % momentum relaxation rate
Gm_pop_em = zeros(1001);
Gm_pop = zeros(1001);

Gpop_abs = zeros(1001); % scattering rate
Gpop_em = zeros(1001);
Gpop_tot = zeros(1001);
wo = (0.03536*e)/hbar; % longitudinal optical phonon energy hbar*wo (eV) (hbar already defined, so explicitly defining wo)
hwo = 0.03536*e; % longitudinal optical phonon energy (eV)

No = 1/(exp(hwo/kT) - 1); 

for i = 1:1001   
   Gpop_abs(i) = real(((q^2)*wo*(ko/kinf - 1))/(2*pi*ko*Eo*hbar*sqrt(2*E(i)/mdos)) * (No*asinh(E(i)/hwo)^(1/2)));
   
   Gm_pop(i) = real(((q^2)*wo*(ko/kinf - 1))/(4*pi*ko*Eo*hbar*sqrt(2*E(i)/mdos)) * (No*sqrt(1 + hwo/E(i)) + (No+1)*sqrt(1 - hwo/E(i)) - hwo*No/E(i) * asinh(E(i)/hwo)^(1/2) + hwo*(No+1)/E(i) * asinh(E(i)/hwo - 1)^(1/2)));
   
   if (E(i) > hbar*wo)
       Gpop_em(i) = ((q^2)*wo*(ko/kinf - 1))/(2*pi*ko*Eo*hbar*sqrt(2*E(i)/mdos)) * ((No + 1)*asinh((E(i)/hwo) -1)^(1/2));
   else
       Gpop_em(i) = 0;
   end
   Gpop_tot(i) = Gpop_abs(i) + Gpop_em(i);
end

figure();
plot(E,Gpop_em(1:1001),'DisplayName','Emission');
hold on;
plot(E,Gpop_abs(1:1001),'DisplayName','Absorption');
plot(E,Gpop_tot(1:1001),'DisplayName','Em + Abs');
title('Polar Optical Phonon Scattering - GaAs');
ylabel('Scattering Rate (1/s)');
xlabel('Energy (J)');
legend;
hold off;

figure();
plot(E,Gpop_tot(1:1001),'DisplayName','\Gamma');
hold on;
plot(E,Gm_pop(1:1001),'DisplayName','\Gamma_m');
legend;
title('POP - Relaxation Rate vs Momentum Relaxation Rate');
ylabel('Scattering Rate (1/s)');
xlabel('Energy (J)');
hold off;

% ---------------------------------------------------------------------------------------
% 3) Intervalley Scattering - absorption and emission for each
% ---------------------------------------------------------------------------------------
% G, 1 valley
% L, 4 valleys (1/2 * 4 sides)
% X, 3 valleys (1/2 * 6 sides)

%Dif_XL = 2e8; % (eV/cm)
Dif_XL = 2e10*e;% (J/m) 
Ein_XL = 0.058*e; % 0.055, 0.041, 0.017 (J)

% Intervalley Deformation Potential
% D_GL = 10e8; % (eV/cm)
% D_GX = 10e8;
% D_LL = 10e8;
% D_LX = 5e8;
% D_XX = 7e8;
D_GL = 10e10*e; % (J/m)
D_GX = 10e10*e;
D_LL = 10e10*e;
D_LX = 5e10*e;
D_XX = 7e10*e;

% Intervalley Phonon Energy
E_GL = 0.0278*e; % (J)
E_GX = 0.0299*e;
E_LL = 0.0290*e; 
E_LX = 0.0293*e;
E_XX = 0.0299*e;

% omegas
w_GL = E_GL/hbar;
w_GX = E_GX/hbar;
w_LL = E_LL/hbar;
w_LX = E_LX/hbar;
w_XX = E_XX/hbar;

% Energy separation between valleys
delta_E_LG = 0.29*e; % (J)
delta_E_XG = 0.45*e; % Gamma to X, but I'm afraid to change it now
delta_E_XL = delta_E_LG - delta_E_XG;
% Number of final valleys
Zf_GL = 4;
Zf_GX = 3;
Zf_LL = 3; % because you can only scatter off of 3
Zf_LX = 3;
Zf_XX = 2; % because you can only scatter off of 2

Zf_LG = 1;
Zf_XG = 1;
Zf_XL = 4;

% Bose-Einstein Dist.
No_GL = 1/(exp(E_GL/kT) - 1);
No_GX = 1/(exp(E_GX/kT) - 1);
No_LL = 1/(exp(E_LL/kT) - 1);
No_LX = 1/(exp(E_LX/kT) - 1);
No_XX = 1/(exp(E_XX/kT) - 1);

% Primitives
gc3d_abs_GL = zeros(1001);
gc3d_em_GL = zeros(1001);
Giv_abs_GL = zeros(1001);
Giv_em_GL = zeros(1001);
Giv_GL = zeros(1001);

gc3d_abs_LG = zeros(1001);
gc3d_em_LG = zeros(1001);
Giv_abs_LG = zeros(1001);
Giv_em_LG = zeros(1001);
Giv_LG = zeros(1001);

gc3d_abs_GX = zeros(1001);
gc3d_em_GX = zeros(1001);
Giv_abs_GX = zeros(1001);
Giv_em_GX = zeros(1001);
Giv_GX = zeros(1001);

gc3d_abs_XG = zeros(1001);
gc3d_em_XG = zeros(1001);
Giv_abs_XG = zeros(1001);
Giv_em_XG = zeros(1001);
Giv_XG = zeros(1001);

gc3d_abs_LL = zeros(1001);
gc3d_em_LL = zeros(1001);
Giv_abs_LL = zeros(1001);
Giv_em_LL = zeros(1001);
Giv_LL = zeros(1001);

gc3d_abs_LX = zeros(1001);
gc3d_em_LX = zeros(1001);
Giv_abs_LX = zeros(1001);
Giv_em_LX = zeros(1001);
Giv_LX = zeros(1001);

gc3d_abs_XL = zeros(1001);
gc3d_em_XL = zeros(1001);
Giv_abs_XL = zeros(1001);
Giv_em_XL = zeros(1001);
Giv_XL = zeros(1001);

gc3d_abs_XX = zeros(1001);
gc3d_em_XX = zeros(1001);
Giv_abs_XX = zeros(1001);
Giv_em_XX = zeros(1001);
Giv_XX = zeros(1001);

for i = 1:1001
   gc3d_abs_GL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GL - delta_E_LG));
   gc3d_em_GL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GL - delta_E_LG));
   
   gc3d_abs_LG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GL + delta_E_LG));
   gc3d_em_LG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GL + delta_E_LG));
   
   gc3d_abs_GX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GX - delta_E_XG));
   gc3d_em_GX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GX - delta_E_XG));
   
   gc3d_abs_XG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GX + delta_E_XG));
   gc3d_em_XG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GX + delta_E_XG));
   
   gc3d_abs_LL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_LL));
   gc3d_em_LL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_LL));
   
   gc3d_abs_LX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_LX - delta_E_XL));
   gc3d_em_LX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_LX - delta_E_XL));
   
   gc3d_abs_XL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_LX + delta_E_XL));
   gc3d_em_XL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_LX + delta_E_XL));
   
   gc3d_abs_XX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_XX));
   gc3d_em_XX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_XX));
   
   
   Giv_abs_GL(i) = (pi*D_GL^2*Zf_GL)/(2*rho*w_GL) * No_GL * gc3d_abs_GL(i); %absorption
   Giv_em_GL(i) = (pi*D_GL^2*Zf_GL)/(2*rho*w_GL) *(No_GL + 1) * gc3d_em_GL(i); %emission
   Giv_GL(i) = Giv_abs_GL(i) + Giv_em_GL(i);
   
   Giv_abs_LG(i) = (pi*D_GL^2*Zf_LG)/(2*rho*w_GL) * No_GL * gc3d_abs_LG(i); %absorption
   Giv_em_LG(i) = (pi*D_GL^2*Zf_LG)/(2*rho*w_GL) *(No_GL + 1) * gc3d_em_LG(i); %emission
   Giv_LG(i) = Giv_abs_LG(i) + Giv_em_LG(i);
   
   Giv_abs_GX(i) = (pi*D_GX^2*Zf_GX)/(2*rho*w_GX) * No_GX * gc3d_abs_GX(i); 
   Giv_em_GX(i) = (pi*D_GX^2*Zf_GX)/(2*rho*w_GX) *(No_GX + 1) * gc3d_em_GX(i); 
   Giv_GX(i) = Giv_abs_GX(i) + Giv_em_GX(i);
   
   Giv_abs_XG(i) = (pi*D_GX^2*Zf_XG)/(2*rho*w_GX) * No_GX * gc3d_abs_XG(i); 
   Giv_em_XG(i) = (pi*D_GX^2*Zf_XG)/(2*rho*w_GX) *(No_GX + 1) * gc3d_em_XG(i); 
   Giv_XG(i) = Giv_abs_XG(i) + Giv_em_XG(i);
   
   Giv_abs_LL(i) = (pi*D_LL^2*Zf_LL)/(2*rho*w_LL) * No_LL * gc3d_abs_LL(i);
   Giv_em_LL(i) = (pi*D_LL^2*Zf_LL)/(2*rho*w_LL) *(No_LL + 1) * gc3d_em_LL(i);
   Giv_LL(i) = Giv_abs_LL(i) + Giv_em_LL(i);
   
   Giv_abs_LX(i) = (pi*D_LX^2*Zf_LX)/(2*rho*w_LX) * No_LX * gc3d_abs_LX(i); 
   Giv_em_LX(i) = (pi*D_LX^2*Zf_LX)/(2*rho*w_LX) *(No_LX + 1) * gc3d_em_LX(i); 
   Giv_LX(i) = Giv_abs_LX(i) + Giv_em_LX(i);
   
   Giv_abs_XL(i) = (pi*D_LX^2*Zf_XL)/(2*rho*w_LX) * No_LX * gc3d_abs_XL(i); 
   Giv_em_XL(i) = (pi*D_LX^2*Zf_XL)/(2*rho*w_LX) *(No_LX + 1) * gc3d_em_XL(i); 
   Giv_XL(i) = Giv_abs_XL(i) + Giv_em_XL(i);
   
   Giv_abs_XX(i) = (pi*D_XX^2*Zf_XX)/(2*rho*w_XX) * No_XX * gc3d_abs_XX(i); 
   Giv_em_XX(i) = (pi*D_XX^2*Zf_XX)/(2*rho*w_XX) *(No_XX + 1) * gc3d_em_XX(i); 
   Giv_XX(i) = Giv_abs_XX(i) + Giv_em_XX(i);
end

figure();
plot(E,Giv_GL(1:1001),'DisplayName','\Gamma-L');
hold on;
plot(E,Giv_LG(1:1001),'DisplayName','L-\Gamma');
title('Intervalley Scattering - GaAs - \Gamma, L');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');
legend;
hold off;

figure();
plot(E,Giv_GX(1:1001),'DisplayName','\Gamma-X');
hold on;
plot(E,Giv_XG(1:1001),'DisplayName','X-\Gamma');
title('Intervalley Scattering - GaAs - \Gamma, X');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');
legend;
hold off;

figure();
plot(E,Giv_LL(1:1001));
title('Intervalley Scattering - GaAs - L-L');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');

figure();
plot(E,Giv_LX(1:1001),'DisplayName','L-X');
hold on;
plot(E,Giv_XL(1:1001),'DisplayName','X-L');
title('Intervalley Scattering - GaAs - L,X');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');
legend;
hold off;

figure();
plot(E,Giv_XX(1:1001));
title('Intervalley Scattering - GaAs - X-X');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');

% ---------------------------------------------------------------------------------------
% 4) Ionized impurity scattering at doping densities of 10^17 and 10^19 cm^-3
% ---------------------------------------------------------------------------------------
 
% Nd_17 = 1e17; % cm^-3
% Nd_19 = 1e19;
Nd_17 = 1e23; % m^-3
Nd_19 = 1e25;

k = zeros(1001);
for i = 1:1001
    k(i) = sqrt(2*mdos*E(i)/hbar^2);
end

Ld_17 = sqrt((Ep*kT)/(e^2 * Nd_17)); % Debye Length
Ld_19 = sqrt((Ep*kT)/(e^2 * Nd_19));

Gi_17 = zeros(1001);
Gi_19 = zeros(1001);
Gmi_17 = zeros(1001); % momentum relaxation rate
Gmi_19 = zeros(1001);
for i = 1:1001
   Gi_17(i) = (Nd_17 * Z^2 * e^4 * Ld_17^4 * mdos)/(hbar^3 * pi * Eo^2 * Es^2) * ((k(i))/(4*(k(i))^2*Ld_17^2 + 1));
   Gi_19(i) = (Nd_19 * Z^2 * e^4 * Ld_19^4 * mdos)/(hbar^3 * pi * Eo^2 * Es^2) * ((k(i))/(4*(k(i))^2*Ld_19^2 + 1));
   
   % Momentum Relaxation Rate
   Gmi_17(i) = (Nd_17*(Z^2)*(e^4)*mdos)/(8*pi*(hbar^3)*(Eo^2)*(Es^2)*(k(i)^3)) * (log(1+4*(k(i)^2)*(Ld_17)^2) - (4*(k(i)^2)*(Ld_17^2))/(1+4*(k(i)^2)*(Ld_17^2)));
   Gmi_19(i) = (Nd_19*(Z^2)*(e^4)*mdos)/(8*pi*(hbar^3)*(Eo^2)*(Es^2)*(k(i)^3)) * (log(1+4*(k(i)^2)*(Ld_19)^2) - (4*(k(i)^2)*(Ld_19^2))/(1+4*(k(i)^2)*(Ld_19^2)));
end

figure();
plot(E(10:1001),Gi_17(10:1001),'DisplayName','Nd = 10^1^7 cm^-^3');
hold on;
plot(E(10:1001),Gi_19(10:1001),'DisplayName','Nd = 10^1^9 cm^-^3');
title('Ionized Impurity Scattering - GaAs');
ylabel('Scattering Rate (1/s)');
xlabel('Energy (J)');
legend;
hold off;

figure();
plot(E(10:1001),Gi_17(10:1001),'DisplayName','\Gamma - Nd = 10^1^7 cm^-^3');
hold on;
title('Ionized Impurity - Relaxation Rate vs. Momentum Relaxation Rate');
ylabel('Scattering Rate (1/s)');
xlabel('Energy (J)');
plot(E(10:1001),Gmi_17(10:1001),'DisplayName','\Gamma_m - Nd = 10^1^7 cm^-^3');
legend;
hold off;

figure();
plot(E,Gi_19(1:1001),'DisplayName','\Gamma - Nd = 10^1^9 cm^-^3');
hold on;
title('Ionized Impurity - Relaxation Rate vs. Momentum Relaxation Rate');
ylabel('Scattering Rate (1/s)');
xlabel('Energy (J)');
plot(E,Gmi_19(1:1001),'DisplayName','\Gamma_m - Nd = 10^1^9 cm^-^3');
legend;
hold off;
% ---------------------------------------------------------------------------------------
% Part B) Relaxation Rates and Momentum Relaxation Rates for POP and
% Ionized Scattering in the same range -- included in respective scattering
% sections
% ---------------------------------------------------------------------------------------



  