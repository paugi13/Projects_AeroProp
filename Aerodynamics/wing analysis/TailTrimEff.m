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
TailEff = 0.721;    % to be changed.
tailIncidence = input('Tail Incidence: ');
flightReg = input('Flight regime? (Takeoff / Cruise) (1/0): ');  % results do not depend on rho and V.

numSt = buildStringAD(tailIncidence);

% Geometry data
XcgMTOW = 14.4408977318597;
Xwing = 15.0317525;
WingTailD = 12.991;
VolumeCoeff = 1.05;
XtailCG = -((Xwing-XcgMTOW)+WingTailD);

AnalysisSize = 5000;
% Flap deflection vector
DE_flap = linspace(-50, 40, AnalysisSize); % flap deflection (deg, positive:down)

% Iteration data
ToL = 500; % tolerated error [Nm]
FuselageAoA = linspace(-15, 20, 100);
tailTrimAngle = zeros(1, length(FuselageAoA));

for i = 1:length(FuselageAoA)
    avAngle = 0;
    for j = 1:length(DE_flap)
        [ResMoment, TailS, CMalpha] = PitchingEffectiveness(XcgMTOW, ...
            WingTailD,Xwing, XtailCG, flightReg, tailIncidence, ...
            tailCLSlope, VolumeCoeff, FuselageAoA(i), TailEff, DE_flap(j));
        if ResMoment < ToL
            avAngle = 1;
            aux = j;
            break
        end
    end
    
    if avAngle == 1
        tailTrimAngle(i) = DE_flap(aux);
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

print(fig1, directSave, '-dpdf', '-r0', '-bestfit');