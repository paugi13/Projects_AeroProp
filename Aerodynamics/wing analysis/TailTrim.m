%% Calculate needed flap angle to trim the wing
% Load tail parameters
clc;
clear;
close all;

loadTailVersion = input('Tail Version to Test?: ');
directLoad = join(['wing analysis/workspaces/TailParameters', ...
    num2str(loadTailVersion)]);
load(directLoad);

% Independent horizontal tail data
A0p = [ 0 0 ]; % root and tip section zero-lift angles (deg)
CM0p = [ 0.0005927 0.0001866 ]; % root and tip section free moments
CDP = [ 0.00594 0.000203 0.00543;% root section CD0, k1 and k2  (airfoil CD curve)
    0.00554 3.75e-5 0.00664 ] ;  % tip section CD0, k1 and k2
DE25 = 20;

% Geometry data
XcgMTOW = 14.4408977318597;
Xwing = 15.0317525;
WingTailD = 13.65;
VolumeCoeff = 1.05;
XtailCG = -((Xwing-XcgMTOW)+WingTailD);

% Flap deflection vector
DE_flap = linspace(-40, 10, 500); % flap deflection (deg, positive:down)
FlapCorr = 1; % flap effectiviness (>=1) because it is seen as a correction

% Iteration data
N = 25 ; % number of panels along the span
ALPHA = [ -10. 20] ; % angles of attack for analysis (deg) 
ToL = 5000; % tolerated error [Nm]
FuselageAoA = linspace(0, 20, 50);
tailTrimAngle = zeros(1, length(FuselageAoA));
flightReg = 1;  % results do not depend on rho and V.

for i = 1:length(FuselageAoA)
    avAngle = 0;    % boolean for available trim angle
    for j = 1:length(DE_flap)
        [cl_local, c4nods, force_coeff, chord, MAC] = GetSolution(N, ...
            ALPHA, FlapCorr, YF_pos, CF_ratio, DE_flap(j), A0p, CM0p, ...
            CDP, AR, TR, DE25, ETIP);
        tailClSlope = polyfit(ALPHA, force_coeff(7, :), 1);
        [ResMoment, TailS, CMalpha] = PitchingCalc(XcgMTOW, WingTailD,...
            Xwing, XtailCG, flightReg, tailIncidence, tailClSlope,...
            VolumeCoeff, FuselageAoA(i));
        if ResMoment < ToL
            avAngle = 1;
            aux = j;
            break
        end
    end
    
    if avAngle == 1
        tailTrimAngle(i) = DE_flap(j);
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

