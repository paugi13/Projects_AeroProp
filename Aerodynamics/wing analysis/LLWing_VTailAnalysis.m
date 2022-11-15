% -------------------------------------------------------------------------     
% VERTICAL TAIL MAIN ANALYSIS WITH DATA SAVING.
% -------------------------------------------------------------------------     
clc; clear; close all;
format long;
addpath(genpath(fileparts(mfilename('fullpath'))));

% -------------------------------------------------------------------------
%% INPUT DATA
% -------------------------------------------------------------------------

% Wing planform (assumes planar wing)

AR = 4;   % aspect ratio
TR = 0.35 ;   % taper ratio (raiz y cola)
DE25 = 20; % sweep angle at c/4 (deg)

ETIP = -0.5; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip)
%-2.245
A0p = [ 0 0 ]; % root and tip section zero-lift angles (deg)
CM0p = [ 0.0005927 0.0001866 ]; % root and tip section free moments
CDP = [ 0.00594 0.000203 0.00543;% root section CD0, k1 and k2  (airfoil CD curve)
    0.00554 3.75e-5 0.00664 ] ;  % tip section CD0, k1 and k2
% Depending on reynolds number

%% FLAP CONFIGURATION
opt = input('Add flap to vertical tail? Y/N (1/0): ');
tailIncidence = 0; % trim the engine at normal flight 
printPlots = input('Print results? Y/N (1/0): ');

switch opt
    case 1
        YF_pos = [ 0.0 1.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
        CF_ratio = 0.35 ;  % flap_chord/chord ratio
        DE_flap = -40; % flap deflection (deg, positive:down)
        FlapCorr = 1; % flap effectiviness (>=1) because it is seen as a correction
    case 0
        YF_pos = [ 0.0 0]; % 2y/b initial and final position of the flap/aileron in the half-wing
        CF_ratio = 0 ;  % flap_chord/chord ratio
        DE_flap = 0; % flap deflection (deg, positive:down)
        FlapCorr = 1.0 ; % flap effectiviness (<=1)
end

%% Simulation data (by the time being only longitudinal analysis)
N = 100 ; % number of panels along the span
ALPHA = [ -10. -9.0 -8.0 -7.0 -5.0 -4.0 -2.0 0. 0.5 0.75 1.0 3.0 5.0 5.25 5.5 6.0 ...
    7.0 9.0 10 10.5  14 20] ; % angles of attack for analysis (deg) 

% -------------------------------------------------------------------------
%% LIFTING LINE SOLUTION
% -------------------------------------------------------------------------
[cl_local, c4nods, force_coeff, chord, MAC] = GetSolution(N, ALPHA, FlapCorr, ...
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
title("\textbf{$Lift contributions$ vs. $Tail span$}");
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
plot(c4nods(2,:), cl_local(:, 5), 'b', 'LineWidth', 1);
plot(c4nods(2,:), cl_local(:, 7), 'r', 'LineWidth', 1);
xlabel("$x/b$ $\left[\mathrm{-}\right]$");
ylabel("$C_l$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
legend('$\alpha = -5^{\circ}$', '$\alpha = -2^{\circ}$', 'location', 'south');
hold off

%% Moment compensation
% Engine: GE Passport 2013. 
rudderCLSlope = polyfit(ALPHA, force_coeff(7,:), 1);

saveRes = input('Save rudder parameters? Y/N (1/0): ');
if saveRes == 1
    confVersion = input('Rudder Version: ');
    directSave = join(['wing analysis/workspaces/RudderParameters', ...
    num2str(confVersion)]);
    save(directSave, 'tailIncidence', 'ETIP', 'YF_pos', 'rudderCLSlope',...
        'AR', 'TR', 'CF_ratio', 'DE_flap');
end

%% Print plots
if printPlots == 1
    if size(ALPHA, 2) > 10
        if CF_ratio == 0
            print(fig1, 'wing analysis/plots/simpleTail_CL_Alpha', ...
                '-dpdf', '-r0', '-bestfit');
            print(fig2, 'wing analysis/plots/simpleTail_CMLE_CL', ...
                '-dpdf', '-r0', '-bestfit');
            print(fig3, 'wing analysis/plots/simpleTail_Basic', ...
                '-dpdf', '-r0', '-bestfit');
            print(fig4, 'wing analysis/plots/simpleTail_LD_alpha', ...
                '-dpdf', '-r0', '-bestfit');
            print(fig5, 'wing analysis/plots/simpleTail_LOCAL_Cl', ...
                '-dpdf', '-r0', '-bestfit');
        else
            print(fig1, 'wing analysis/plots/FlapTail_CL_Alpha', ... 
                '-dpdf', '-r0', '-bestfit');
            print(fig2, 'wing analysis/plots/FlapTail_CMLE_CL', ...
                '-dpdf', '-r0','-bestfit');
            print(fig3, 'wing analysis/plots/FlapTail_Basic', ...
                '-dpdf', '-r0', '-bestfit');
            print(fig4, 'wing analysis/plots/FlapTail_LD_alpha', ...
                '-dpdf', '-r0', '-bestfit');
        end
    end
end