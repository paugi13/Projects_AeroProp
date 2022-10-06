%% Propulsion project: Part 1
clear all;
%% Known parameters
% Engine: Spey RSp.4 Mk.511-8

% Geomètriques
Weng = 1126*9.81;
eixos = 2;
n_lpc = 5;
n_hpc = 12;
n = n_hpc+n_lpc;
alpha = 0.64;
T_cr=50790;

% Performance
SFC_cr = 1.7e-5;

%% REVISAR
m_cr = 0;
OPR_0 = 18.4;
Tt4 = 1350;
%1350
%% -----------

%Del problema
h_0 = 0;
M_cr = 0;
M0 = M_cr;
m_ssl=89.4;

%ISA
P0 = 101325;
T0 = 288.15;
g = 9.81;
%% -------------------

pi_s = OPR_0^(1/n);
pi_lpc = pi_s^(n_lpc);
pi_hpc = pi_s^(n_hpc);

% Físiques 
R = 287;
gamma_c = 1.4;
gamma_t = 1.33;
cpt = gamma_t/(gamma_t-1)*R;
cpc = gamma_c/(gamma_c-1)*R;
H = 43000000;

theta0 = 1+ (gamma_c-1)/2*M0^2;
delta0 = theta0^(gamma_c/(gamma_c-1));

Tt0 = T0 * theta0;

% Dades donades problema i no pel fabricant
pi_d = 0.98;
rend_fan = 0.89;
rend_lpc = 0.89; %Fan
rend_hpc = 0.89; % 
pi_b = 0.96;
rend_b = 0.99;
rend_hpt = 0.91;
rend_lpt = 0.93;
rend_mH = 0.993;
rend_mL = 0.997;
pi_np = 0.99;
pi_ns = 0.99;

% In this case we can assume pi_fan = pi_lpc because there isn't LPC.
pi_fan = pi_lpc;


%% Necessary equations

tau_lpc = 1 + (pi_lpc^((gamma_c-1)/gamma_c)-1)/rend_lpc;
tau_hpc = 1 + (pi_hpc^((gamma_c-1)/gamma_c)-1)/rend_hpc;
tau_fan = 1 + (pi_fan^((gamma_c-1)/gamma_c)-1)/rend_fan;

Tt3 = T0 * theta0 * tau_lpc * tau_hpc;
f = (cpt*Tt4 - cpc*Tt3) / (rend_b*H);

Tt2 = Tt0;
tau_hpt = 1-(1/rend_mH)*(1/(1+f))*(cpc/cpt)*(Tt2/Tt4)*tau_lpc*(tau_hpc-1);
tau_lpt = 1-(1/rend_mL)*(1/(1+f))*(cpc/cpt)*(Tt2/Tt4)*(1/tau_hpt)*((tau_lpc-1)+alpha*(tau_fan-1));

Tt9 = Tt4 * tau_lpt * tau_hpt;

pi_hpt = (1-((1-tau_hpt)/rend_hpt))^(gamma_t/(gamma_t-1));
pi_lpt = (1-((1-tau_lpt)/rend_lpt))^(gamma_t/(gamma_t-1));

%Pt19/P0
rel_press_fan = delta0*pi_d*pi_fan*pi_ns;
Pt19 = P0*rel_press_fan;

M19options = [1 sqrt(2/(gamma_c-1)*(rel_press_fan^((gamma_c-1)/gamma_c)-1))];
M19=min(M19options);
%Pt9/P0
rel_press_p9 = delta0*pi_d*pi_lpc*pi_hpc*pi_b*pi_hpt*pi_lpt*pi_np;
M9options = [1 sqrt(2/(gamma_t-1)*(rel_press_p9^((gamma_t-1)/gamma_t)-1))];
M9=min(M9options);

Pt9=rel_press_p9*P0;

%càlculs previs pel thrust
T9=Tt9/(1+M9^2*(gamma_t-1)/2);
T19=Tt0*tau_fan/(1+M19^2*(gamma_c-1)/2);


a0=sqrt(gamma_c*R*T0);
u9=M9*sqrt(gamma_t*R*T9);
u19=M19*sqrt(gamma_c*R*T19);

p9=Pt9*(2/(gamma_t+1))^(gamma_t/(gamma_t-1));
p19=Pt19*(2/(gamma_c+1))^(gamma_c/(gamma_c-1));

%Specific thrust.
Thrust_esp=((1+f)*u9+(1+f)*(R*T9)/(u9)*(1-P0/p9)-M0*a0) + (alpha*(u19-M0*a0)+alpha*(R*T19)/(u19)*(1-P0/p19));

% Coreflow in %
coreflow=((1+f)*u9+(1+f)*(R*T9)/(u9)*(1-P0/p9)-M0*a0)/Thrust_esp*100;
bypass=(alpha*(u19-M0*a0)+alpha*(R*T19)/(u19))*(1-P0/p19)/Thrust_esp*100;

% Specific consumption
Cs=f/Thrust_esp;
Is = 1/(Cs*g);

m=T_cr/Thrust_esp;
mf=m*f;
m_bypass=m*alpha;
m_tot=m_bypass+m;


