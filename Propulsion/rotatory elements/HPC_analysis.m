close all
clear
clc

rh = 0.32;
rt = 0.48;
rm = (rh+rt)/2;
Mi = 0.4;
AT_max = 45;
alpha_1 = 30;
Tti = 329.7244;
P1 = 27446.35;
Pt1 = P1*1.4802*2.3551;
Ti = Tti/(1+(1.4-1)/2*Mi^2);
pi_hpc = 7.8129;
tau_hpc = 1.898;
gamma_c = 1.4;

%% 1

mu_p = (gamma_c-1)/gamma_c*log(pi_hpc)/log(tau_hpc);
mu_p_r = (1+mu_p)/2;

%% 2

Tt1= Ti*(1+(gamma_c-1)/2*Mi^2);
Tt_diff = Tt1*(tau_hpc-1);
etapes_min = Tt_diff/AT_max;
et_min = 7;
AT_real = Tt_diff/et_min;

v1 = Mi*sqrt(gamma_c*Ti*287)*sind(alpha_1);
Cp = 287*gamma_c/(gamma_c-1);

syms x
eq = Cp*AT_real == x^2 - x*2*v1;
U = solve(eq,x);
U = double(U(2));


%% 3
% 5th stage.
% V1_5 is only at the middle section and it isnt constant along the radius

M1_5 = 0.4*sqrt(Ti/(Ti+4*AT_real));
T1_5 = Ti + 4*AT_real;
Tt1_5 = Tt1 + 4*AT_real;

% Necessary calculus to obtain hub and tip radius at 5th stage. 
% rm stays the same (repeated stages), but the values of rh and rt can't be
% assumed the same as 1st stage. 
Pt1_5 = Pt1 * (Tt1_5/Tt1)^(gamma_c*mu_p/(gamma_c-1));
P1_5 = P1*Pt1_5/Pt1*((1+(gamma_c-1)/2*Mi^2)/((1+(gamma_c-1)/2*M1_5^2)))...
    ^(gamma_c/(gamma_c-1));

dens_1 = P1/(287*Ti);
dens_5 = P1_5/(287*T1_5);
h1 = (rt-rm)*2;
h5 = h1*dens_1/dens_5;

rt = rm + h5/2;
rh = rm - h5/2;

T2_1 = Ti + (1-0.5)*AT_real;
T2_5 = Ti + (5-0.5)*AT_real;
omega = U/rm;
Uh = omega*rh;
Ut = omega*rt;
V1_5 = M1_5*sqrt(gamma_c*287*T1_5);


% rt5 = rm + h5/2;
% rh5 = rm - h5/2;
% Using radial distributions theory
v1_5_m = V1_5*sind(alpha_1);
k1 = v1_5_m*rm;
u1_5_m = V1_5*cosd(alpha_1);

% Root
v1_5_h = k1/rh;
% alpha_1_h = asind(v1_5_h/V1_5); 
% u1_5_h = V1_5*cosd(alpha_1_h);
beta_1_h = atand((Uh-v1_5_h)/u1_5_m);

% Tip
v1_5_t = k1/rt;
% alpha_1_t = asind(v1_5_t/V1_5);
% u1_5_t = V1_5*cosd(alpha_1_t);
beta_1_t = atand((Ut-v1_5_t)/u1_5_m);

% ct = U*(U-aux*sind(30));



