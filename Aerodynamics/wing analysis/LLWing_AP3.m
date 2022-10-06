
% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodinàmica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

clc; clear all; close all;
format long;

% -------------------------------------------------------------------------
% INPUT DATA
% -------------------------------------------------------------------------

% Wing planform (assumes planar wing)

AR = 21.3 ;   % aspect ratio
TR = 1/5.55 ;   % taper ratio (raiz y cola)
DE25 = 15.7 ; % sweep angle at c/4 (deg)

ETIP = -3.8; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip)

A0p = [ -1 -1 ]; % root and tip section zero-lift angles (deg)
CM0p = [ 0.012 0.012 ]; % root and tip section free moments
CDP = [ 0.0066 -0.0061 0.0063;% root section CD0, k1 and k2  (airfoil CD curve)
    0.01 -0.0075 0.0096 ] ;  % tip section CD0, k1 and k2

% Flap/aileron (symmetrical deflection)

YF_pos = [ 0.0 0.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2 ;  % flap_chord/chord ratio
DE_flap = 0.0; % flap deflection (deg, positive:down)
FlapCorr = 1.0 ; % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)

N = 100 ; % number of panels along the span

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10.] ; % angles of attack for analysis (deg) 

% -------------------------------------------------------------------------
% LIFTING LINE SOLUTION
% -------------------------------------------------------------------------

% Wing discretization (lenghts are dimensionless with the wing span)

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

% Assembly of the influence coefficients matrix (needed once)

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

% Solve circulations for different angles of attack

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

% Loads calculation using plane Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

% -------------------------------------------------------------------------
% POSTPROCESS ( your work starts here...)
% -------------------------------------------------------------------------


%% Part 1
%Datos
%momento de cabeceo 
pi=3.14;
k=0.75;
bf=20;
alpha=[-10*pi/180,20*pi/180];
CMf=258.2996*alpha/(3.02*0.5*(100/3.6)^2*1.225);
plot(alpha,CMf)
title 'Coeficiente de momento de cabeceo'
xlabel '\alpha[rad]'
ylabel 'M_F [Nm]'
%constantes
a=0.3;
b=0.6;
lh=2;
SB=pi*(b/2)^2;
Ss=3.06;
Re=3.91e6;
Cf=0.074/(Re^(1/5));
q=0.5*1.225*(100/3.6)^2;
Vf=3.02;
% ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10.];
S_bal = pi*0.3.*sqrt((1+(tand(ALPHA)).^2).*(1./(1/0.09+(tand(ALPHA)).^2)));
figure
plot(ALPHA,S_bal);

%Drag parasito
%dfb=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss/SB);
%db=sqrt(Sb/0.7854);
%Cdb=0.029*(db/b)^3/sqrt(Cdfb);
Cd0=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss/SB);


%% Part 2. 
%1. OK
%drag parasito en función de alfa
Cd01=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss./S_bal);

%momento de cabeceo en función de alfa
Mf=2*k*q*Vf*ALPHA*pi/180;
C_mf=Mf/(q*Vf);
C_mf2 = (pi*k*q*2*0.09*ALPHA*pi/180*4/3)/(q*Vf);
figure
plot(ALPHA,C_mf)
title 'Momento de cabeceo'
xlabel '\alpha [º]'
ylabel 'C_{M_f}'
grid on

%2.Wing's CL. OK
figure
plot(ALPHA, force_coeff(7,:))
xlabel('\alpha [º]')
ylabel('C_L')
grid on
%alpa0=-4.87
%Clalpha=0.08

%3.Wing's CM_LE - CL
CM_le=force_coeff(5,:);
figure
plot(force_coeff(7,:),CM_le, 'b');
xlabel('C_L')
ylabel('C_{M_{LE}}')
grid on
cm_0 = (CM_le(3)-CM_le(2))/(force_coeff(7,3)-force_coeff(7,2))*(-force_coeff(7,2))+CM_le(2);
x_ca_c = (CM_le(3)-CM_le(2))/(force_coeff(7,3)-force_coeff(7,2));
X_ac = -1*x_ca_c*mac*bf;
 
%4. Additional / Basic Lift. OK
Cla=(cl_local(:,2)-cl_local(:,3))/(force_coeff(7,2)-force_coeff(7,3));
Clb=cl_local(:,3)-Cla*force_coeff(7,3);
figure
plot(Cla,'b')
hold on
plot(Clb,'r')
xlabel('Paneles')
ylabel('Cl')
legend('Cla','Clb')
grid on
hold off

%5. Cdo - k model. OK
v = linspace(-1.5, 1, 150);
pol = polyval([0.0211 -0.0063613 0.10948], v);
CBW= Cd01+force_coeff(11,:);
figure
plot(force_coeff(7,:),CBW)
hold on
plot(v, pol, 'r');
grid on
xlabel('Cl')
ylabel('Cd_{BW}')
hold off

%6. CMcg-CL
area_trap = 18.8; 
X_cg=1.3;
Cm_cg=cm_0-force_coeff(7,:)*(X_ac - X_cg)/(mac*bf)+C_mf2; 
figure
plot(force_coeff(7,:),Cm_cg);
title 'Coeficiente de momento respecto el centro de gravedad'
xlabel 'Cl'
ylabel 'Cm_{cg}'
axis equal
grid on 


%% Part 3. Analysis. 
%Datos
S_wing = 18.8;
S_total = S_wing + Ss;
weight = 9.81*20*S_total;
dens = 1.225;

%1. Stability margin
mean_aero = mac*bf;
marg = 100*(X_ac - X_cg)/mean_aero;

%2. Trim coindition 
%sweep = 15.7  (if the margin is calculated well)
%washout = -3.8

%3. CL_max 
%Usando distribuciones de lift básico y adicional 

%Compute Cl_max (local) as a Reynolds function using XFLR5
chord_real = chord*bf;
%Re de cada sección
Re_section = 1.225*100*chord_real/(3.6*1.789e-5);
%Hay una relación exponencial entre Re y cl_local maximo. 
cl_max_vec = log(Re_section/1e-6)/16.58;
CL_max_vec = (cl_max_vec - Clb)./Cla;
CL_min = min(CL_max_vec);

%Stall speed. Calculating lift equaling weight with CL_min. 
u_stall = sqrt(2*weight/(dens*S_wing*CL_min));


%4. CD0 - K Model
%Done in ap2











