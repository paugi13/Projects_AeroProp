%% Calculate needed flap angle to trim the wing
% Load tail parameters

clc;
clear;
close all;

% Load tail without flap parameters. 
load('wing analysis/workspaces/TailParametersPlainTail');

% Tail Parameters
tailCLSlope = polyfit(ALPHA, CLTail, 1);
TailEff = 0.6;    % to be changed.
tailIncidence = input('Tail Incidence: ');

% Geometry data
XcgMTOW = 14.4408977318597;
Xwing = 15.0317525;
WingTailD = 13.65;
VolumeCoeff = 1.05;
XtailCG = -((Xwing-XcgMTOW)+WingTailD);

% Flap deflection vector
DE_flap = linspace(-40, 10, 500); % flap deflection (deg, positive:down)

% Iteration data
ToL = 5000; % tolerated error [Nm]
FuselageAoA = linspace(0, 20, 50);
tailTrimAngle = zeros(1, length(FuselageAoA));
flightReg = 1;  % results do not depend on rho and V.

for i = 1:length(FuselageAoA)
    avAngle = 0;
    for j = 1:length(DE_flap)
        [ResMoment, TailS, CMalpha] = PitchingEffectiveness(XcgMTOW, ...
            WingTailD,Xwing, XtailCG, flightReg, tailIncidence, ...
            tailCLSlope, VolumeCoeff, FuselageAoA, TailEff);
        if ResMoment < ToL
            avAngle = 1;
            aux = j;
        end
    end
    
    if avAngle == 1
        tailTrimAngle(i) = DE_flap(aux);
    else
        error('Unable to find trim angle for %s.', FuselageAoA(i));
    end
end

%% POSTPROCESS
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
directSave = join(['wing analysis/plots/TailTrim', ...
    num2str(loadTailVersion)]);

fig1 = figure(1);
hold on
title("\textbf{$\delta_e$ vs $\alpha_{wb}$ $i_w = 5.25^\circ$}");
plot(FuselageAoA, tailTrimAngle, 'b', 'LineWidth', 1);
xlabel("$\alpha_{wb}$ $\left[\mathrm{^\circ}\right]$");
ylabel("$\delta_e$ $\left[\mathrm{^\circ}\right]$");
grid on;
grid minor;
box on;
hold off

print(fig1, directSave, '-dpdf', '-r0', '-bestfit');