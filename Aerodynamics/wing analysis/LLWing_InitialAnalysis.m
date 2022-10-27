
% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodinàmica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

clc; clear; close all;
format long;
addpath(genpath(fileparts(mfilename('fullpath'))));

% -------------------------------------------------------------------------
%% INPUT DATA
% -------------------------------------------------------------------------

% Wing planform (assumes planar wing)

AR = 8.6 ;   % aspect ratio
TR = 0.35 ;   % taper ratio (raiz y cola)
DE25 = 22.5; % sweep angle at c/4 (deg)

ETIP = -4; % tip twist (deg, negative for washout)

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
ALPHA = [ -10. -9.0 -8.0 -7.0 -5.0 -4.0 -2.0 0. 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10 10.5 11 12 14 20] ; % angles of attack for analysis (deg) 

% -------------------------------------------------------------------------
%% LIFTING LINE SOLUTION
% -------------------------------------------------------------------------
[cl_local, c4nods, force_coeff] = GetSolution(N, ALPHA, FlapCorr, ...
    YF_pos, CF_ratio, DE_flap, A0p, CM0p, CDP, AR, TR, DE25, ETIP);

% -------------------------------------------------------------------------
%% POSTPROCESS
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

eff = force_coeff(7,:)./force_coeff(11,:);
Cla=(cl_local(:,2)-cl_local(:,3))/(force_coeff(7,2)-force_coeff(7,3));
Clb=cl_local(:,3)-Cla*force_coeff(7,3);
CM_le=force_coeff(5,:); 

% Wing's Cl - alpha
fig1 = figure(1);
hold on
title("\textbf{$C_L$ vs. $\alpha$}");
% title("\textbf{Plot $C_L$ vs. $\alpha$ $b_f/b = 2/3$}");
plot(ALPHA, force_coeff(7,:), 'b', 'LineWidth', 1)
xlabel("$\alpha$ $\left[\mathrm{^\circ}\right]$");
ylabel("$C_L$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
hold off

disp(force_coeff(7, end));

% Wing's CM_LE - CL. OK
fig2 = figure(2);
hold on
title("\textbf{$C_{M_{LE}}$ vs. $\alpha$}");
plot(force_coeff(7,:),CM_le, 'b');
xlabel("$C_L$ $\left[\mathrm{-}\right]$");
ylabel("$C_{M_{LE}}$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
hold off

% Additional / Basic Lift. OK
fig3 = figure(3);
hold on
title("\textbf{$Lift contributions$ vs. $Wing span$}");
plot(Cla, 'b', 'LineWidth', 1)
plot(Clb, 'r', 'LineWidth', 1)
xlabel("Spanwise station $\left[\mathrm{-}\right]$");
ylabel("$C_{l}$ $\left[\mathrm{-}\right]$");
legend('$C_{la}$','$C_{lb}$')
grid on;
grid minor;
box on;
hold off

% Efficiency (Fuselage drag is missing)
fig4 = figure(4);
hold on
title("\textbf{$L/D$ vs. $\alpha$}");
plot(ALPHA, eff, 'b', 'LineWidth', 1)
xlabel("$\alpha$ $\left[\mathrm{^\circ}\right]$");
ylabel("$L/D$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
hold off

fig5 = figure(5);
hold on
title("\textbf{Local $C_l$ vs. Spanwise station }");
plot(c4nods(2,:), cl_local(:, 12), 'b', 'LineWidth', 1);
plot(c4nods(2,:), cl_local(:, 17), 'r', 'LineWidth', 1);
xlabel("$x/b$ $\left[\mathrm{-}\right]$");
ylabel("$C_l$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
legend('$\alpha = 5^{\circ}$', '$\alpha = 10^{\circ}$', 'location', 'south');
hold off

if size(ALPHA, 2) > 10
    if CF_ratio == 0
        print(fig1, 'wing analysis/plots/simpleWing_CL_Alpha', '-dpdf', '-r0', ...
            '-bestfit');
        print(fig2, 'wing analysis/plots/simpleWing_CMLE_CL', '-dpdf', '-r0', ...
            '-bestfit');
        print(fig3, 'wing analysis/plots/simpleWing_Basic', '-dpdf', '-r0', ...
            '-bestfit');
        print(fig4, 'wing analysis/plots/simpleWing_LD_alpha', '-dpdf', '-r0', ...
            '-bestfit');
        print(fig5, 'wing analysis/plots/simpleWing_LOCAL_Cl', '-dpdf', '-r0', ...
            '-bestfit');
    else
        print(fig1, 'wing analysis/plots/FlapWing_CL_Alpha', '-dpdf', '-r0', ...
            '-bestfit');
        print(fig2, 'wing analysis/plots/FlapWing_CMLE_CL', '-dpdf', '-r0', ...
            '-bestfit');
        print(fig3, 'wing analysis/plots/FlapWing_Basic', '-dpdf', '-r0', ...
            '-bestfit');
        print(fig4, 'wing analysis/plots/FlapWing_LD_alpha', '-dpdf', '-r0', ...
            '-bestfit');
    end
end









