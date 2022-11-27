% Code for fuselage aerodynamic performance analysis.
clc;
clear;
close all;

WantedAoA = input('Wanted cruise AoA? ');

% Load cruise configuration
load('wing analysis/workspaces/wingLiftdist525');

%% Flight data
% ALPHA = [ -10. -8.0 -6 -4.0 -2 0. 1 2 3 4.0 8.0 10.] ;
ALPHA = -10:0.5:10;

for i = 1:length(ALPHA)
    if ALPHA(i) == WantedAoA
        pos = i;
        a = 1;
        break
    end
end

if a~=1
    error('Angle of attack not in the alpha vector');
end

%%  Drag analysis
% Fuselage dimensions.
r = 1.589;  % fuselage radius
b = 2*r;    % fuselage diameter
lh = 31.26320;  % fuselage length

Re=10e6;
Cf=0.074/(Re^(1/5));

Ss = 219.62;    %   Total wetted area (ellipsoid area)
S_bal = pi*r.*sqrt((1+(tand(ALPHA)).^2).*(1./(1/r^2+(tand(ALPHA)).^2)));
% Frontal area as a function of AoA.

Cd0 = Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss./S_bal);

%% Pitching moment coefficient
Vf = 133.76;
v = 857/3.6;
rho = 0.363918;
q = 0.5*rho*v^2;
k = 0.9;

% function 4*((3.178/2)^2(1-(x^2/(31.2632/2)^2))) 
% integration limits -lh/2 -> lh/2
intResult = 52.62474*4;
Mf = pi/2*k*q*ALPHA*pi/180*intResult;
Cmf = Mf/(q*Vf);

%% Drag and pitching moment computation. 
wingS = 65.258;
MfCruise = Cmf(pos)*q*Vf;
Df = Cd0(pos)*q*(Vf^(2/3));
Dw = wingCD*q*wingS;

DTotal = Df + Dw;
LTotal = wingCL*q*wingS;

E = LTotal/DTotal;


%% POTSPROCESS
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig1 = figure(1);
hold on
title("\textbf{$C_{D_0}$ vs. $\alpha$}");
plot(ALPHA, Cd0, 'b', 'LineWidth', 1)
xlabel("$\alpha$ $\left[\mathrm{^\circ}\right]$");
ylabel("$C_{D_0}$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
hold off

fig2 = figure(2);
hold on
title("\textbf{$C_{M_f}$ vs. $\alpha$}");
plot(ALPHA, Cmf, 'b', 'LineWidth', 1)
xlabel("$\alpha$ $\left[\mathrm{^\circ}\right]$");
ylabel("$C_{M_f}$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
hold off

print(fig1, 'wing analysis/plots/FuselageCD0', '-dpdf', '-r0', '-bestfit');
print(fig2, 'wing analysis/plots/FuselageCMF', '-dpdf', '-r0', '-bestfit');