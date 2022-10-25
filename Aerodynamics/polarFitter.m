clc;
clear; 
close all;

% Script for polar reading and fitting.
addpath(genpath(fileparts(mfilename('fullpath'))));

% genData = table2array(readtable('SC20414-Cl_Cd.csv'));
genData = readmatrix('SC20410-Cl_alpha.csv');

% For Re(1.000.000) and Re(3.000.000), 1-2, 5-6.
fig1= figure(1);
hold on
title("\textbf{Plot $C_l$ vs. $\alpha$}");
plot(genData(1:end,1), genData(1:end,2), 'b' , 'LineWidth', 1);
plot(genData(1:end,7), genData(1:end,8), 'r' , 'LineWidth', 1);
xlabel("$\alpha$ $\left[\mathrm{^\circ}\right]$");
ylabel("$C_l$ $\left[\mathrm{-}\right]$");
legend('Re = 1.000.000', 'Re = 3.000.000', 'location', 'northwest');
grid on;
grid minor;
box on;
hold off

print(fig1, 'airfoil design/plots/SC(2)-0410_alpha_Cl', '-dpdf', '-r0', ...
    '-bestfit');


