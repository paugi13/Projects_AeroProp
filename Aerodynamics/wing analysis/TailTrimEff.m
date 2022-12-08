%% Plane pitching moment trim 
% Calculate needed flap angle to trim the wing
clc;
clear;
close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

% Load tail without flap parameters. 
load('wing analysis/workspaces/TailParametersPlainTail');

% Tail Parameters
tailCLSlope = polyfit(ALPHA, CLTail, 1);
TailEff = 0.67;    % xh/c = 0.7
tailIncidence = input('Tail Incidence: ');
flightReg = input('Flight regime? (Takeoff / Cruise) (1/0): ');  % results do not depend on rho and V.

numSt = buildStringAD(tailIncidence);

% Geometry data
XcgMTOW = 13.3101;  % XcgOEW = 14.2753   XcgMTOW = 13.3101
Xwing = 13.9085;        % 13.9085;
WingTailD = 15.01; 
VolumeCoeff = 1.15; %1.05;
XtailCG = -((Xwing-XcgMTOW)+WingTailD);

AnalysisSize = 10000;
% Flap deflection vector
DE_flap = linspace(-50, 40, AnalysisSize); % flap deflection (deg, positive:down)

% Iteration data
ToL = 500; % tolerated resultant moment [Nm]
FuselageAoA = linspace(-15, 20, 500);
tailTrimAngle = zeros(1, length(FuselageAoA));
M1V = zeros(1, length(FuselageAoA));
M2V = zeros(1, length(FuselageAoA));
M3V = zeros(1, length(FuselageAoA));
epsV = zeros(1, length(FuselageAoA));

auxVec = zeros(1, length(FuselageAoA));
for i = 1:length(FuselageAoA)
    avAngle = 0;
    for j = 1:length(DE_flap)
        [M1, M2, M3, ResMoment, TailS, CMalpha, eps] = PitchingEffectiveness(XcgMTOW, ...
            WingTailD,Xwing, XtailCG, flightReg, tailIncidence, ...
            tailCLSlope, VolumeCoeff, FuselageAoA(i), TailEff, DE_flap(j));
        if abs(ResMoment) < ToL
            avAngle = 1;
            aux = j;
            break
        end
    end
    
    if avAngle == 1
        tailTrimAngle(i) = DE_flap(aux);
        auxVec(i) = ResMoment;
        M1V(i) = M1;
        M2V(i) = M2;
        M3V(i) = M3;
        epsV(i) = eps;
    else
        error('Unable to find trim angle for %s.', FuselageAoA(i));
    end
end

disp(tailTrimAngle(1));
%% POSTPROCESS
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

if flightReg == 1
directSave = 'wing analysis/plots/TailTrimFlap';
else
    directSave = 'wing analysis/plots/TailTrimCruise';
end

fig1 = figure(1);
hold on
title("\textbf{$\delta_e$ vs $\alpha_{wb}$ $i_{hw} = " + tailIncidence ...
    + "^\circ$}");
plot(FuselageAoA, tailTrimAngle, 'b', 'LineWidth', 1);
xlabel("$\alpha_{wb}$ $\left[\mathrm{^\circ}\right]$");
ylabel("$\delta_e$ $\left[\mathrm{^\circ}\right]$");
grid on;
grid minor;
box on;
hold off
% 2.12, 8.48 
print(fig1, directSave, '-dpdf', '-r0', '-bestfit');


% XcgMTOW = 14.5612;  % XcgOEW = 14.2753   XcgMTOW = 13.8263
% Xwing = 15.1585;
% WingTailD = 13.76;