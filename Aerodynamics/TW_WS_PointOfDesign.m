clear;
clc;
close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%dades generals
S = 76.17 ; %superficie [m^2]
mtow = 39848.43 ; %kg
v_cruise = 842.0*1000/3600 ; %velocitat de creuer [m/s]
T_to = 125e3 ; % takeoff thrust [N]
AR = 8.94 ; %aspect ratio [-]
g = 9.81;


%% Similar aircrafts' point of design

a220H = 0.347498475;
a220V = 5309.542386;
crj900H = 0.337218527;
crj900V = 5325.269792;
crj700H = 0.367366812;
crj700V = 4864.816181;

 
%% Takeoff
% takeoff at??
sigma = 0.974/1.225 ;
cLmax = 1.8 ;
% Takeoff parameter
TOP = mtow*9.81/S * mtow*9.81/T_to * 1/cLmax * 1/sigma ; 

%% Cruise
rho = 0.363918 ; %kg/m^3; 11km
q = 0.5*rho*v_cruise^2 ; %dynamic pressure 
W_to = mtow*9.81 ;
W_cr = W_to * 0.85 ; %la relació ha de ser aprox. 0.85 segons els apunts que tenen a greva
e = 0.8 ; %oswald efficiency factor
T_cr = 22241.108 ; %N (T en creuer és molt diferent que en takeoff, dada trobada a internet -web guardada a l'excel-)

cD = T_cr / (q*S) ; %T = D
cL = W_cr / (q*S) ; %L = W
cD_0 = cD - cL^2 / (pi*e*AR) ;


%% Distància d'aterratge
Vstall = 76.56 ; %m/s (calculada a partir de les mitges de les Vstall calculades de cada aeornau semblant)
W_S_1 = 0.5*Vstall^2*0.974*2.5 ; %línia vertical del plot (hem agafat 126m/s com a referència perq per tema flaps, etc no serà 76.56)

%% 2n segment 
gamma = atan(1.2/100) ;
N_e = 2 ;
W_2 = W_to ;
T_2seg = T_to ; 
q_2 = 0.5*1.21*(1.2*Vstall)^2 ;
cD_2 = (T_2seg/2 + W_2*gamma) / (q_2*S) ; %El thrust aquí és la meitat ja que hi ha fallada crítica d'un dels motors
cL_2 = W_2 / (q_2*S) ;
%T_W_2seg = N_e/(N_e - 1) * T_to/2*1/(T_2seg/2) * W_2/W_to * (cD_2/cL_2 +
%gamma) --> utilitzem la fórmula d'abaix perq la relació de pesos i thrust
%és 1 (l'equació utilitzada està treta de: https://www.fzt.haw-hamburg.de/pers/Scholz/dimensionierung/2seg_nn.htm
T_W_2seg = N_e/(N_e - 1) * (cD_2/cL_2 + gamma) ;

%% càlculs dels vectors per plotejar
W_S = linspace(100,8500,100) ;
T_W_to = zeros(length(W_S),1) ;
T_W_cr = zeros(length(W_S),1) ;
for nElem = 1:length(W_S)
    T_W_to(nElem) = W_S(nElem) / TOP * 1/cLmax * 1/sigma ;
    T_W_cr(nElem) = T_to/T_cr * W_cr/W_to * q/(W_S(nElem)*W_cr/W_to) * ( cD_0 + (W_S(nElem)*W_cr/W_to)^2 * 1/(q^2*pi*AR*e)) ;
end

%% Point of Design
finalMTOW = 39850*g;
CLMax = 3.062;
CLPoD = 2.1;
PointOfDesignVer = 0.5*Vstall^2*0.974*CLPoD;
PointOfDesignHor = 0.374;

finalT = finalMTOW*PointOfDesignHor;
finalS = finalMTOW/PointOfDesignVer;


%% Diagrama
figure
hold on
plot(W_S,T_W_to,W_S,T_W_cr)
xline(W_S_1, 'r');
yline(T_W_2seg)
% xline(PointOfDesignVer, 'b');
% yline(PointOfDesignHor, 'b');

scatter(PointOfDesignVer, PointOfDesignHor, 'r', 'filled');
scatter(a220V, a220H,   'filled');
scatter(crj700V, crj700H, 'filled');
scatter(crj900V, crj900H,'filled');
plot([0 PointOfDesignVer], [PointOfDesignHor PointOfDesignHor], '--', 'color', 'r' );
plot([PointOfDesignVer PointOfDesignVer], [0 PointOfDesignHor],  '--', 'color', 'r');
xlim([0 8500])
ylim([0 1])
legend('Takeoff', 'Cruise', 'Landing distance', '2nd segment climb', 'Design point', 'A220', 'CRJ700', 'CRJ900', 'location', 'north')
title('T/W - W/S')
xlabel('$\frac{W}{S} [N/m^2]$', 'interpreter','latex')
ylabel('$\frac{T}{W} [-]$','interpreter','latex')
grid on
hold off
