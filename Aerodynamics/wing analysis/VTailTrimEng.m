%% Plane yaw moment trim. 
% Failing engine criteria. 
clc; clear; close all;
format long;
addpath(genpath(fileparts(mfilename('fullpath'))));

load('wing analysis/workspaces/RudderParameters0');
T = 78930;  % get exact value 78930 72250;
% yFusRoot = 1.374;
% yEngRoot = 2.8;   % to be defined
yTotal = -3.7534;  

% wing characteristics
wingS = 65.258;
b = 2*11.85;
% rudder geometry
lv = 12.2;  % distance 
XcgMTOW = 14.4408977318597;
Xwing = 15.0317525;
Vv = 0.11;
tailS = Vv*b*wingS/lv;
eff = 0.81;
BvBr = 1;

% flight profile: 80% of stall speed
v = 196.777/3.6*1.13;   % § 25.1513	Minimum control speed.
rho = 1.225;
q = 0.5*rho*v^2;
mu = 0.97;  % relation between alpha 

AnalysisSize = 500;
TVector = linspace(0, T, AnalysisSize);
TrimAngle = zeros(1, length(TVector));
ToL = 10;

for i = 1:length(TVector)
    avAngle = 0;
    TMoment = TVector(i)*yTotal;
    CnCoeff = rudderCLSlope(1)*Vv*mu*eff*BvBr;
    TrimAngle(i) = -TMoment/(CnCoeff*wingS*b*q);
end

disp(['Tail Area = ', num2str(tailS)]);

%% POSTPROCESS

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% directSave = join(['wing analysis/plots/TailTrim', ...
%     numSt]);

fig1 = figure(1);
hold on
title("\textbf{$\delta_e$ vs $\alpha_{wb}$ $y_{CG}= " + yTotal + "~m$}");
plot(TVector, TrimAngle, 'b', 'LineWidth', 1);
xlabel("$T$ $\left[\mathrm{N}\right]$");
ylabel("$\delta_r$ $\left[\mathrm{^\circ}\right]$");
grid on;
grid minor;
box on;
xlim([0 T]);
hold off

print(fig1, 'wing analysis/plots/VTailTrim', '-dpdf', '-r0',...
    '-bestfit');
