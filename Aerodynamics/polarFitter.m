clc;
clear; 
close all;

% Script for polar reading and fitting.
addpath(genpath(fileparts(mfilename('fullpath'))));

% genData = table2array(readtable('SC20414-Cl_Cd.csv'));
genData = readmatrix('SC20414-Cl_Cd.csv');

% For Re(1.000.000) and Re(3.000.000), 1-2, 5-6.
figure
hold on
title("\textbf{Plot $C_d$ vs. $C_l$}");
plot(genData(1:end-20,1), genData(1:end-20,2), 'b' , 'LineWidth', 1);
plot(genData(1:end-20,5), genData(1:end-20,6), 'r' , 'LineWidth', 1);
xlabel("$C_l$ $\left[\mathrm{-}\right]$");
ylabel("$C_d$ $\left[\mathrm{-}\right]$");
legend('Re = 1.000.000', 'Re = 3.000.000', 'location', 'northwest');
grid on;
grid minor;
box on;
hold off


