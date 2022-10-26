clc;
clear; 
close all;

% Script for polar reading and fitting.
addpath(genpath(fileparts(mfilename('fullpath'))));

% genData = table2array(readtable('SC20414-Cl_Cd.csv'));
genData = readmatrix('SC20410-Cl_Cd.csv');

% For Re(1.000.000) and Re(3.000.000), 1-2, 5-6.
fig1= figure(1);
hold on
title("\textbf{Plot $C_d$ vs. $C_l$}");
plot(genData(20:end-50,1), genData(20:end-50,2), 'b' , 'LineWidth', 1);
plot(genData(20:end-50,7), genData(20:end-50,8), 'r' , 'LineWidth', 1);
xlabel("$C_l$ $\left[\mathrm{^\circ}\right]$");
ylabel("$C_d$ $\left[\mathrm{-}\right]$");
legend('Re = 1.000.000', 'Re = 3.000.000', 'location', 'northwest');
grid on;
grid minor;
box on;
hold off

print(fig1, 'airfoil design/plots/SC(2)-0410_Cl_Cd', '-dpdf', '-r0', ...
    '-bestfit');


