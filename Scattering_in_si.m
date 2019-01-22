% Si

close all;
clear;

e = 1.602e-19; % charge of an electron/proton
%hbar = 6.582e-16; % Reduced Planck's Constant (eV*s)
hbar = 1.054e-34; % (J*s)
%rho = 2.329; % mass density (g/cm^3)
rho = 2329; % (kg/m^3)
%kT = 0.0259; % (eV) 
kT = 0.0259*e; %(J)
%vs = 9.04e5; % longitudinal acoustic velocity (cm/s)
vs = 9040; % (m/s)
%Dac = 9.5; % electron acoustic deformation potential (eV)
Dac = 9.5 * e; % (J)
%mo = 9.11e-28; % g 
mo = 9.11e-31; % kg
ml = 0.98*mo;
mt = 0.19*mo;
%mdos = 6^(2/3)*(ml*mt^2)^(1/3);
mdos = (0.98*0.19^2)^(1/3)*mo;
Ec = 0;
E = linspace(0,2*e,1001); % 0 to 2*e Joules

% ---------------------------------------------------------------------------------------
% Acoustic Phonon Scattering - within elastic and equipartition approximations
% ---------------------------------------------------------------------------------------

gc3d = zeros(1001);
G = zeros(1001);

for i = 1:1001
   %gc3d(i) = 6/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * sqrt(E(i) - Ec);
   %G(i) = 2 * 2*pi/hbar * Dac^2 * kT/(2*rho*vs^2) * (1/2) * gc3d(i);
   G(i) = (Dac^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*mdos/(hbar^2))^(3/2) * E(i)^(1/2);
end

plot(E,G(1:1001));
title('Acoustic Phonon Scattering - Si');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');

% ---------------------------------------------------------------------------------------
% Equivalent Intervalley Scattering - g and f processes, absorption and emission for each
% ---------------------------------------------------------------------------------------

% f processes - 1 final valley in Si
% g processes - 4 final valleys in Si
Giv_abs_g = zeros(1001);
Giv_em_g = zeros(1001);
Giv_both_g = zeros(1001);

Giv_abs_f = zeros(1001);
Giv_em_f = zeros(1001);
Giv_both_f = zeros(1001);

gc3d_abs_g = zeros(1001);
gc3d_em_g = zeros(1001);
gc3d_abs_f = zeros(1001);
gc3d_em_f = zeros(1001);

% g-type X-X
%Dif_TA_g = 0.5e8; %(eV/cm)
%Dif_TA_g = 5e9; % (eV/m)
Dif_TA_g = 5e9 * e; % (J/m)
%Dif_LA_g = 0.8e8;
%Dif_LA_g = 8e9; % (eV/m)
Dif_LA_g = 8e9 * e; % (J/m)
%Dif_LO_g = 11e8;
%Dif_LO_g = 110e9; %(eV/m)
Dif_LO_g = 110e9 * e; %(J/m)
Dif_g = (Dif_TA_g + Dif_LA_g + Dif_LO_g) / 3;

% f-type X-X
%Dif_TA_f = 0.3e8;
%Dif_TA_f = 3e9; %(eV/m)
Dif_TA_f = 3e9 * e; %(J/m)
%Dif_LA_f = 2e8;
%Dif_LA_f = 20e9; %(eV/m)
Dif_LA_f = 20e9 * e; %(J/m)
%Dif_LO_f = 2e8;
%Dif_LO_f = 20e9; %(eV/m)
Dif_LO_f = 20e9 * e; %(J/m)
Dif_f = (Dif_TA_f + Dif_LA_f + Dif_LO_f) / 3;

Zf_f = 1; % Number of final valleys
Zf_g = 4;

delta_Eif = 0; % 0 for equivalent valley scattering in Si
Eiv_g = e*(0.012 + 0.019 + 0.062) / 3; % (J) %0.062? 
Eiv_f = e*(0.019 + 0.047 + 0.059) / 3;

wif_f = Eiv_f/hbar;
wif_g = Eiv_g/hbar;

Nif_f = 1/(exp(Eiv_f/kT) - 1); % given by Bose-Einstein Distribution Function
Nif_g = 1/(exp(Eiv_g/kT) - 1);

for i = 1:1001
   
   gc3d_abs_f(i) = 6/(2*pi^2) * (2*ml/hbar^2)^(3/2) * real(sqrt(E(i) + Eiv_f - delta_Eif));
   %gc3d_em_f(i) = 6/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_f - delta_Eif));
   
   gc3d_abs_g(i) = 6/(2*pi^2) * (2*mt/hbar^2)^(3/2) * real(sqrt(E(i) + Eiv_g - delta_Eif));
   %gc3d_em_g(i) = 6/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_g - delta_Eif));
   
   Giv_abs_f(i) = (pi*Dif_f^2*Zf_f)/(2*rho*wif_f) * Nif_f * gc3d_abs_f(i); %absorption
   %Giv_em_f(i) = (pi*Dif_f^2*Zf_f)/(2*rho*wif_f) *(Nif_f + 1) * gc3d_em_f(i); %emission
   
   Giv_abs_g(i) = (pi*Dif_g^2*Zf_g)/(2*rho*wif_g) * Nif_g * gc3d_abs_g(i);
   %Giv_em_g(i) = (pi*Dif_g^2*Zf_g)/(2*rho*wif_g) *(Nif_g + 1) * gc3d_em_g(i);
   
   if (E(i) >= Eiv_f)
       gc3d_em_f(i) = 6/(2*pi^2) * (2*ml/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_f - delta_Eif));
       Giv_em_f(i) = (pi*Dif_f^2*Zf_f)/(2*rho*wif_f) *(Nif_f + 1) * gc3d_em_f(i); %emission
   else
       Giv_em_f(i) = 0;
   end
   
   if (E(i) >= Eiv_g)
       gc3d_em_g(i) = 6/(2*pi^2) * (2*mt/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_g - delta_Eif));
       Giv_em_g(i) = (pi*Dif_g^2*Zf_g)/(2*rho*wif_g) *(Nif_g + 1) * gc3d_em_g(i);
   else
       Giv_em_g(i) = 0;
   end
   
   Giv_both_f(i) = Giv_abs_f(i) + Giv_em_f(i);
   Giv_both_g(i) = Giv_abs_g(i) + Giv_em_g(i);
   
end

figure();
plot(E,Giv_both_f(1:1001));
title('Intervalley Scattering - f-process');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');

figure();
plot(E,Giv_both_g(1:1001));
title('Intervalley Scattering - g-process');
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');

% ---------------------------------------------------------------------------------------
% Ionized impurity scattering at doping densities of 10^17 and 10^19 cm^-3
% ---------------------------------------------------------------------------------------
q = e;
%Nd_17 = 1e17; % cm^-3
Nd_17 = 1e23; % m^-3
%Nd_19 = 1e19;
Nd_19 = 1e25; % m^-3
Z = 1; 
Es = 11.7; % relative dielectric constant
Eo = 8.854e-12; %F/m vacuum permittivity
Ep = Es*Eo; 
a = 5.43e-10; %m
%k = linspace(0,sqrt(2*mdos*2/hbar^2),1001);
k = zeros(1001);
for i = 1:1001
    k(i) = sqrt(2*mdos*E(i)/hbar^2);
end

Ld_17 = sqrt((Ep*kT)/(q^2 * Nd_17)); % Debye Length
Ld_19 = sqrt((Ep*kT)/(q^2 * Nd_19));

Gi_17 = zeros(1001);
Gi_19 = zeros(1001);

for i = 1:1001
   Gi_17(i) = (Nd_17 * Z^2 * q^4 * (Ld_17^4) * mdos)/(hbar^3 * pi * Eo^2 * Es^2) * (k(i)/(4*(k(i)^2)*(Ld_17^2) + 1));
   Gi_19(i) = (Nd_19 * Z^2 * q^4 * (Ld_19^4) * mdos)/(hbar^3 * pi * Eo^2 * Es^2) * (k(i)/(4*(k(i)^2)*(Ld_19^2) + 1));
end

figure();
plot(E(10:1001),Gi_17(10:1001),'DisplayName','Nd = 10^1^7 cm^-^3');
hold on;
plot(E(10:1001),Gi_19(10:1001),'DisplayName','Nd = 10^1^9 cm^-^3');
hold off;
legend;
xlabel('Energy (J)');
ylabel('Scattering Rate (1/s)');
title('Ionized Impurity Scattering - Si');