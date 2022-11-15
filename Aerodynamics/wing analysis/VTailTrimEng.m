%% Plane yaw moment trim. 
% Failing engine criteria. 

load('wing analysis/workspaces/RudderParameters0');
T = 80000;  % get exact value 
yEng = -5;   % to be defined

% wing characteristics
wingS = 65.258;
b = 2*11.85;
% rudder geometry
lv = 12.5;  % distance 
XcgMTOW = 14.4408977318597;
Xwing = 15.0317525;
Vv = 0.1;
Sv = Vv*b*wingS/lv;
maxDef = 25;
eff = 0.4;

engMoment = T*yEng;

% flight profile: 80% of stall speed
v = 196.777/3.6*0.8;
rho = 1.225;
q = 0.5*rho*v^2;
mu = 0.97;  % relation between alpha 

AnalysisSize = 500;
AngleSize = 500;
TVector = linspace(0, T, AnalysisSize);
rudderAngleVector = linspace(0, maxDef, AngleSize);
TrimAngle = zeros(1, length(TVector));
ToL = 10;

for i = 1:length(TVector)
    avAngle = 0;
    TMoment = TVector(i)*yEng;
    CnCoeff = rudderCLSlope(1)*Vv*mu*eff;
    syms angle
    eq = -TMoment == CnCoeff*wingS*b*angle;
    sol = solve(eq, angle);
    TrimAngle(i) = double(sol);
end

%% POSTPROCESS

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
directSave = join(['wing analysis/plots/TailTrim', ...
    numSt]);

fig1 = figure(1);
hold on
title("\textbf{$\delta_e$ vs $\alpha_{wb}$ $i_{hw} = -8.42^\circ$}");
plot(TVector*1000, TrimAngle, 'b', 'LineWidth', 1);
xlabel("$T$ $\left[\mathrm{kN}\right]$");
ylabel("$\delta_r$ $\left[\mathrm{^\circ}\right]$");
grid on;
grid minor;
box on;
hold off
