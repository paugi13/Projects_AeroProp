clc;
clear; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% AILERON EFFECTIVENESS ANALYSIS. 
% Code to get aileron effectiveness. 
incr_alpha = 0.5;
n_alpha = -5:incr_alpha:20;
alpha = -5;                                      %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.

%% Define new camber line for supercrital airfoils
coord = table2array(readtable('NASASC(2)0012.csv'));
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

%% Aileron parameters. Flap approach
def = 1;
xhVector = linspace(0, 1, 100);
effVector = zeros(1, length(xhVector));
% NASASC(2)0410 zero lift angle = -2.2917
al0Simple = 0;
al0Vector = zeros(1, length(xhVector));

%% Airfoil analysis
for j = 1:length(xhVector)
    for i=1:length(n_alpha)
        [cl_dist(1, i), cm_dist(1, i)] = cl_cm0_NACA_alpha_function(point,...
            pan, alpha, pos, xhVector(j), def);
        cm_0(1, i) = cl_dist(1, i)*0.25 + cm_dist(1, i);
        if int32(alpha) == 5
           m = (cl_dist(1,i)-cl_dist(1,(i-1)))/(alpha-(alpha-incr_alpha));
           al0Vector(j) = -cl_dist(1,(i-1))/m + (alpha-incr_alpha);   
        end
        alpha = alpha+incr_alpha;
    end
    effVector(j) = abs(al0Simple - al0Vector(j))/def;
    alpha = -5;
end

%% Postprocess
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig1 = figure(1);
hold on
title("\textbf{Elevator effectiveness vs. Hinge position}");
plot(1-xhVector, effVector, 'b', 'LineWidth', 1);
xlabel("$c_a/c$ $\left[\mathrm{-}\right]$");
ylabel("$\tau$ $\left[\mathrm{-}\right]$");
grid on
grid minor
axis equal
xlim([0 1]);
ylim([0 effVector(1)]);
hold off

print(fig1, 'plots/tailEffectiveness', '-dpdf', '-r0', '-bestfit');
