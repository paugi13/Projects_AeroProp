%% Lift computation
clc; clear; close all;
format long;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Load workspace for desired incidence angle
% Only one analysis at the time

flightReg = input('Takeoff or cruise regime? (1/0): ' );
wantedAoA = input('Wing incidence angle?: ');

if mod(wantedAoA, 1) ~= 0
    if mod(wantedAoA*10, 1) ~= 0    % two decimals
        st1 = num2str(fix(wantedAoA));
        st2 = num2str(fix(mod(wantedAoA,1)*10));
        st3 = num2str(mod(mod(wantedAoA,1)*10, 1)*10);
        numSt = join([st1, st2, st3]);
    else % one decimal
        st1 = num2str(fix(wantedAoA));
        st2 = num2str(mod(wantedAoA, 1)*10);
        numSt = join([st1, st2]);
    end
else        % no remainder
    numSt = num2str(wantedAoA);
end

%% Input data
% cruise regime: MTOW?
MTOW = 39850*9.81;

%% Calculations
if flightReg == 1   %takeoff
    TOW = MTOW;
    TOSpeed = 270.392/3.6;
    rho = 0.974;
    q = 0.5*rho*TOSpeed^2;
    direct = join(['wing analysis/workspaces/wingLiftdistFlap', ...
        numSt]);
    load(direct);
else    %cruise
    cruiseWeight = MTOW * 0.85;
    cruiseSpeed = 857/3.6;
    rho = 0.363918;
    q = 0.5*rho*cruiseSpeed^2;
    direct = join(['wing analysis/workspaces/wingLiftdist', ...
        numSt]);
    load(direct);
end

liftDist = 2*q*Cl_Values;
d = 5;      % polynomial degree = d - 1.
dist = @(x) 0;
for i = 1:length(polinomialFit)
    dist = @(x) dist(x) + polinomialFit(i)*x.^d;
    d = d - 1;
end
startCoord = spanCoords(1);
endCoord = spanCoords(end);
TotalLift = 2*q*integral(dist, startCoord, endCoord);

if flightReg == 1
    disp(['TOW =  ', num2str(TOW)]);
    disp(['Incidence angle = ', num2str(wantedAoA)]);
    disp(['Provided Lift = ', num2str(TotalLift)]);
else
    disp(['Cruise weight =  ', num2str(cruiseWeight)]);
    disp(['Incidence angle = ', num2str(wantedAoA)]);
    disp(['Provided Lift = ', num2str(TotalLift)]);
end

%% Plot lift distribution
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
legendStr = "$\alpha = " + wantedAoA + "^{\circ}$"; 

fig1 = figure(1);
hold on
title("\textbf{Local $C_l$ vs. Spanwise station }");
plot(spanCoords, Cl_Values, 'b', 'LineWidth', 1);
xlabel("$x/b$ $\left[\mathrm{-}\right]$");
ylabel("$C_l$ $\left[\mathrm{-}\right]$");
grid on;
grid minor;
box on;
legend(legendStr, 'location', 'south');
hold off

%% Save workspace
if flightReg == 1   % takeoff
direct = join(['wing analysis/workspaces/TakeOffData', numSt]);
save(direct , 'polinomialFit', 'q', 'startCoord', 'endCoord');
else    % cruise
direct = join(['wing analysis/workspaces/CruiseData', numSt]);
save(direct , 'polinomialFit', 'q', 'startCoord', 'endCoord');
end