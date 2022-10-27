clc;
clear;
close all;

AR = 8.6 ;   % aspect ratio
TR = 0.35 ;   % taper ratio (raiz y cola)
DE25 = 22.5; % sweep angle at c/4 (deg)

ETIP = [-2 -2.5 -3 -3.5 -4 -4.5 -5 -5.5 -6]; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip)
%-2.245
A0p = [ -2.9531 -2.245 ]; % root and tip section zero-lift angles (deg)
CM0p = [ -0.09686 -0.07562 ]; % root and tip section free moments
CDP = [ 0.006096 -0.002796 0.006573;% root section CD0, k1 and k2  (airfoil CD curve)
    0.012 -0.00159 0.00573 ] ;  % tip section CD0, k1 and k2
% Depending on reynolds number

%% FLAP CONFIGURATION
opt = input('Add flap to wing? Y/N (1/0): ');

switch opt
    case 1
        YF_pos = [ 0.0 0.67]; % 2y/b initial and final position of the flap/aileron in the half-wing
        CF_ratio = 0.3 ;  % flap_chord/chord ratio
        DE_flap = 35; % flap deflection (deg, positive:down)
        FlapCorr = 1; % flap effectiviness (>=1) because it is seen as a correction
    case 0
        YF_pos = [ 0.0 0]; % 2y/b initial and final position of the flap/aileron in the half-wing
        CF_ratio = 0 ;  % flap_chord/chord ratio
        DE_flap = 0; % flap deflection (deg, positive:down)
        FlapCorr = 1.0 ; % flap effectiviness (<=1)
end

%% Simulation data (by the time being only longitudinal analysis)
N = 100 ; % number of panels along the span
ALPHA = 5.0 ; % angles of attack for analysis (deg) 

clLocalSeries = zeros(size(ETIP, 2), N);
maxAcc = zeros(size(ETIP, 2), 1);
maxPos = zeros(size(ETIP, 2), 1);

for i = 1:size(ETIP, 2)
    [clLocalSeries(i, :), c4nods, force_coeff] = GetSolution(N, ALPHA, FlapCorr, ...
        YF_pos, CF_ratio, DE_flap, A0p, CM0p, CDP, AR, TR, DE25, ETIP(i));
    [maxAcc(i, 1), maxPos(i,1)] =  max(clLocalSeries(i, :));
    if maxPos(i,1) > N/2
        maxPos(i,1) = N - maxPos(i,1);
    end
end

% POSTPROCESS
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

legendStrings = " $\epsilon $ = " + string(ETIP) + "$^\circ$";

fig1 = figure(1);
hold on
title("\textbf{Local $C_l$ vs. Spanwise station }");
for i = 1:size(ETIP, 2)
plot(c4nods(2,:), clLocalSeries(i,:), 'LineWidth', 1);
end
plot(c4nods(2,maxPos), maxAcc(:,1), 'r', 'LineWidth', 1);
plot(c4nods(2,100-maxPos), maxAcc(:,1), 'r', 'LineWidth', 1);
xlabel("$x/b$ $\left[\mathrm{-}\right]$");
ylabel("$C_l$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
legend(legendStrings, 'location', 'south');
hold off

print(fig1, 'wing analysis/plots/ETIP_Comparison', '-dpdf', '-r0', ...
            '-bestfit');

