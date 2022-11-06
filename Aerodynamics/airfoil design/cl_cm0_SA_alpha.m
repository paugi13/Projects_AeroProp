clc;
clear; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% SINGLE AIRFOIL ANALYSIS READING .CSV FILE

incr_alpha = 0.01;
n_alpha = -5:incr_alpha:20;
alpha =-5;                                      %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.

%% Define new camber line for supercrital airfoils
coord = table2array(readtable('NASASC(2)0410.csv'));
point = size(coord,1)/2;
pan = point-1;
pos = zeros(point, 2);

for i=1:size(coord,1)/2
pos(i,1) = coord(i,1);
pos(i,2) = (coord(i,2)+coord(i+point, 2))*0.5;
end

cl_dist = zeros(1, length(n_alpha));
cm_dist = zeros(1, length(n_alpha));
cm_0 = zeros(1, length(n_alpha));

%% Airfoil analysis
for i=1:length(n_alpha)
    [cl_dist(1, i), cm_dist(1, i)] = cl_cm0_NACA_alpha_function(point, pan, alpha, pos, 0, 0);
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

fig1 = figure(1);
hold on
title("\textbf{NASASC(2)-0414}");
plot(pos(:,1), pos(:,2), 'r');
plot(coord(1:point, 1), coord(1:point, 2), 'b');
plot(coord(point+1:point*2, 1), coord(point+1:point*2, 2), 'b');
xlabel("$x/c$ $\left[\mathrm{-}\right]$");
ylabel("$y$ $\left[\mathrm{-}\right]$");
grid on
grid minor
axis equal
hold off

fig2 = figure(2);
hold on
title("\textbf{$C_l$ vs $\alpha$}");
plot(n_alpha, cl_dist(1,:), 'b', 'LineWidth', 1);
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_l$ $\left[\mathrm{-}\right]$')
grid on
grid minor
hold off

fig3 = figure(3);
hold on
title("\textbf{$C_{m_0}$ vs $\alpha$}");
plot(n_alpha, cm_0(1,:), 'b' , 'LineWidth', 1);
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_{m_0}$ $\left[\mathrm{-}\right]$')
grid on
grid minor
axis equal
hold off

fig4 = figure(4);
hold on
title("\textbf{$C_m$ vs $\alpha$}");
plot(n_alpha, cm_dist(1,:), 'b' , 'LineWidth', 1);
xlabel('$\alpha$ $\left[\mathrm{^\circ}\right]$')
ylabel('$C_m$ $\left[\mathrm{-}\right]$')
grid on
grid minor
axis equal
hold off

print(fig1, 'plots/NASASC(2)-0414', '-dpdf', '-r0', '-bestfit');
print(fig2, 'plots/SC(2)-0414_Cl_alpha', '-dpdf', '-r0', '-bestfit');