
% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodinàmica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

clc; clear; close all;
format long;
addpath(genpath('LLWING'));
% -------------------------------------------------------------------------
%% INPUT DATA
% -------------------------------------------------------------------------

% Wing planform (assumes planar wing)

AR = 21.3 ;   % aspect ratio
TR = 1/5.55 ;   % taper ratio (raiz y cola)
DE25 = 17 ; % sweep angle at c/4 (deg)

ETIP = -7.1; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip)

A0p = [ -1 -1 ]; % root and tip section zero-lift angles (deg)
CM0p = [ 0.012 0.012 ]; % root and tip section free moments
CDP = [ 0.0066 -0.0061 0.0063;% root section CD0, k1 and k2  (airfoil CD curve)
    0.01 -0.0075 0.0096 ] ;  % tip section CD0, k1 and k2
% Depending on reynolds number

% Flap/aileron (symmetrical deflection)

YF_pos = [ 0.0 0.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2 ;  % flap_chord/chord ratio
DE_flap = 0.0; % flap deflection (deg, positive:down)
FlapCorr = 1.0 ; % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)

N = 100 ; % number of panels along the span

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ] ; % angles of attack for analysis (deg) 

% -------------------------------------------------------------------------
%% LIFTING LINE SOLUTION
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
%momento de cabeceo 
pi=3.14;
k=0.75;
bf=20;
alpha=[-10*pi/180,20*pi/180];
CMf=258.2996*alpha/(3.02*0.5*(100/3.6)^2*1.225);

figure
hold on
plot(alpha,CMf);
title 'Coeficiente de momento de cabeceo'
xlabel '\alpha[rad]'
ylabel 'M_F [Nm]'
hold off

%constantes
a=0.3;
b=0.6;
lh=2;
SB=pi*(b/2)^2;
% Sb=0.25*SB;
Ss=3.06;
Re=3.91e6;
Cf=0.074/(Re^(1/5));
q=0.5*1.225*(100/3.6)^2;
Vf=3.02;

%Area de b en función de al angulo de ataque.
%S_balpha=pi.*0.3*sqrt((1./(1+(tand(ALPHA).^2)/0.09)).*(1+tand(ALPHA).^2));
%figure
%plot(ALPHA,S_balpha);

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ];
S_bal = pi*0.3.*sqrt((1+(tand(ALPHA)).^2).*(1./(1/0.09+(tand(ALPHA)).^2)));
figure
hold on
plot(ALPHA,S_bal);
hold off
%drag parasito

%dfb=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss/SB);
%db=sqrt(Sb/0.7854);
%Cdb=0.029*(db/b)^3/sqrt(Cdfb);

% Fuselage drag
% Cd0=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss/SB);


%% Part 2. 
%1. OK
%drag parasito en función de alpha
Cd01=Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss./S_bal);

%momento de cabeceo en función de alfa
Mf=2*k*q*Vf*ALPHA*pi/180;
C_mf=Mf/(q*Vf);             %Munk Method -> Da momentos demasiado altos
C_mf2 = (pi*k*q*2*0.09*ALPHA*pi/180*4/3)/(q*Vf);    %Por integración del pdf 3.3. 
figure
plot(ALPHA,C_mf)
title('Momento de cabeceo')
xlabel('\alpha [º]')
ylabel('Cm_f')
grid on

%2.Wing's Cl. OK
figure
plot(ALPHA, force_coeff(7,:))
xlabel('\alpha [º]')
ylabel('Cl')
grid on
%alpa0=-4.87
%Clalpha=0.08

%3.Wing's CM_LE - CL. OK
CM_le=force_coeff(5,:); 
figure
plot(force_coeff(7,:),CM_le, 'b');
xlabel('C_L')
ylabel('Cm_{LE}')
grid on
cm_0 = 0.1197;  %Por interpolación
x_ca_c = -1.498;    %Pendiente de la recta
X_ac = -1*x_ca_c*mac*bf;    %Cálculo a partir de la mean aerodynamic chord.


%NECESSARI?
%-----------------------------------------------
%pendent=0.085
% figure
% plot(ALPHA, force_coeff(5,:))
% xlabel('\alpha [º]')
% ylabel('Cm0')
% grid on
% dCmdCl=-0.817;
% Xac=-dCmdCl*1.55;
%-------------------------------------------


sweep_back= 17; %deg
lamb=0.28/1.55;
%X_A no es la posición del centro aerodinámico. 
X_A= tand(sweep_back)*(bf/6)*((1+2*lamb)/(1+lamb));
 
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
pol = polyval([0.021021 -0.0064986 0.1103], v);
CBW= Cd01+force_coeff(11,:);
figure
plot(force_coeff(7,:),CBW)
hold on
plot(v, pol, 'r');
grid on
xlabel('C_L')
ylabel('CD_{BW}')
hold off

%6. CMcg-CL. OK
area_trap = 18.8; 
% Cm0=0.166;
Xcg=1.3;
%Distancias adimensionalizadas.
Cm_cg=cm_0-force_coeff(7,:)*(X_ac - Xcg)/(mac*bf)+C_mf2;  
figure
plot(force_coeff(7,:),Cm_cg);
title 'Coeficiente de momento respecto el centro de gravedad'
xlabel 'C_L'
ylabel 'CM_{cg}'
axis equal
grid on 

%% Part 3. Analysis. 






