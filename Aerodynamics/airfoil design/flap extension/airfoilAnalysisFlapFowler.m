clc;
clear; close all;
% addpath(genpath('C:\Users\Usuari\Desktop\Q7\Projectes\Aerodinàmica - Propulsió\Projects_AeroProp\Aerodynamics\airfoil design'));
addpath(genpath(fileparts(mfilename('fullpath'))));

incr_alpha = 0.01;
n_alpha = -5:incr_alpha:20;
alpha =-5;                                      %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.

% According to NASA's data with GA(W)-1
alphaStall = 10.5;

%% Define new camber line for supercrital airfoils
% GA(W)-1
coord = table2array(readtable('NASASC(2)0414.csv'));
coord_flap = table2array(readtable('fowler_coordinates.csv'));
point = size(coord,1)/2;
pointFlap = size(coord_flap,1)/2;
pan = point - 1;
pos = zeros(point, 2);
posFlapAux = zeros(pointFlap, 2);
panFlap = pointFlap - 1;
xRel = 1;   % flap starts at the end of the airfoil

for i=1:size(coord,1)/2
pos(i,1) = coord(i,1);
pos(i,2) = (coord(i,2)+coord(i+point, 2))*0.5;
end

for i=1:size(coord_flap,1)/2
posFlapAux(i,1) = coord_flap(i,1);
posFlapAux(i,2) = (coord_flap(i,2)+coord_flap(i+pointFlap, 2))*0.5;
end

cl_dist = zeros(1, length(n_alpha));
cm_dist = zeros(1, length(n_alpha));
cm_0 = zeros(1, length(n_alpha));

%% Create flap fowler camber line with rotation matrix
rotAngle = 35;
posFlap = flapCalculator(rotAngle, posFlapAux, pointFlap, xRel);

%% Rotate flap surfaces for plotting
coordRot = flapCalculator(rotAngle, coord_flap, pointFlap*2, xRel);

%% Build general position vector
panTotal = pan + panFlap;
pointTotal = panTotal + 1;
posGen = totalCamberBuilder(pos,posFlap);

%% Airfoil analysis
for i=1:length(n_alpha)
    point=pan+1;
    [cl_dist(1, i), cm_dist(1, i)] = cl_cm0_NACA_alpha_function(pointTotal, panTotal, alpha, posGen, 0, 0);
    cm_0(1, i) = cl_dist(1, i)*0.25 + cm_dist(1, i);
    if int32(alpha) == 5
       m = (cl_dist(1,i)-cl_dist(1,(i-1)))/(alpha-(alpha-incr_alpha));
       al_0 = -cl_dist(1,(i-1))/m + (alpha-incr_alpha);   
    end
    alpha = alpha+incr_alpha;
end


%% Postprocess
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Plot airfoil with camber line
nClStall = 1551;
ClStallData = 3.625;

figure
hold on
title("\textbf{NASASC(2)-0414} - $35^\circ$ \textbf{Flap Fowler}" );
plot(pos(:,1), pos(:,2), 'r');
plot(coord(1:point, 1), coord(1:point, 2), 'b');
plot(coord(point+1:point*2, 1), coord(point+1:point*2, 2), 'b');
% plot(posFlapAux(:,1)+xRel, posFlapAux(:,2), 'r');
% plot(coord_flap(1:pointFlap, 1)+xRel, coord_flap(1:pointFlap, 2), 'b');
% plot(coord_flap(pointFlap+1:pointFlap*2, 1)+xRel, ...
%     coord_flap(pointFlap+1:pointFlap*2, 2), 'b');
plot(posFlap(:,1), posFlap(:,2), 'r');
plot(coordRot(1:pointFlap, 1), coordRot(1:pointFlap, 2), 'b');
plot(coordRot(pointFlap+1:pointFlap*2, 1), ...
    coordRot(pointFlap+1:pointFlap*2, 2), 'b');
grid on
grid minor
axis equal
xlabel("$x/c$ $\left[\mathrm{-}\right]$");
ylabel("$y$ $\left[\mathrm{-}\right]$");
hold off


figure
hold on
title("\textbf{$C_l$ vs $\alpha$}");
scatter(alphaStall, cl_dist(nClStall), [], [1 0 0], 'filled');
scatter(alphaStall, ClStallData, [], [0 1 0], 'filled');
plot(n_alpha, cl_dist(1,:), 'b', 'LineWidth', 1);
% plot([alphaStall alphaStall], [0 cl_dist(nClStall)], '--', 'color', 'r');
% plot([0 alphaStall], [cl_dist(nClStall) cl_dist(nClStall)], '--', 'color', 'r');
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_l$ $\left[\mathrm{-}\right]$')
grid on
grid minor
legend('Calculated stall point', 'Data stall point' ,'location', 'northwest');
hold off

figure
hold on
title("\textbf{$C_{m_0}$ vs $\alpha$}");
plot(n_alpha, cm_0(1,:), 'b' , 'LineWidth', 1);
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_{m_0}$ $\left[\mathrm{-}\right]$')
grid on
grid minor
axis equal
hold off