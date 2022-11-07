clc;
clear; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

% Code to analyze aileron performance given with its geometry als data.
%% Load workspace
wantedAoA = input('Wing incidence angle?: ');
numSt = buildStringAD(wantedAoA);
direct = join(['wing analysis/workspaces/wingLiftdist', ...
    numSt]);
load(direct);

%% Aileron rolling moment calculations
% Ixx using 737 weight distribution data and considering MTOW.
% %Weight wing = 33.53%
MWing = 13361.705;
MFuselage = 26488.295;
R = 1.589;
b = 11.85*2;
Cr = 4.08;
tAverage = Cr*0.14*0.5;
Ixx = MFuselage*R^2 + 1/12*MWing*(b^2 + tAverage^2);

% Wing characteristics
wingS = 65.258;
lambda = 0.35;
chordF = @(y) Cr*(1 + 2*(lambda - 1)/b*y); 
Ct = chordF(b/2);
upAngle = 30;
downAngle = 25;
meanAngle = (abs(upAngle) + abs(downAngle))/2;

% Aileron characteristics
xh = 0.3;
eff = 0.6566;
innerEdge = 0.7;
outerEdge = 0.95;

% Flight profile data
turnSpeed = 350/3.6;
rho = 0.974;
q = 0.5*rho*turnSpeed^2;

%% Rolling moment integration 
yi = innerEdge*b/2;
yo = outerEdge*b/2;
chordFInt = @(y) y^2/2 + 2/3*(lambda - 1)/b*y^3;
intResult = chordFInt(yo) - chordFInt(yi);

%% Aileron lift sections
innerY = innerEdge*b/2;
outerY = outerEdge*b/2;

% Build integral function.
d = 5;      % polynomial degree = d - 1.
dist = @(x) 0;
for i = 1:length(polinomialFit)
    dist = @(x) dist(x) + polinomialFit(i)*x.^d;
    d = d - 1;
end

resultInt = integral(dist, innerY, outerY);

% Surface integral
surfFunc = polyfit(spanCoords, wingChord(51:end,1)', 5);

d = 5;      % polynomial degree = d - 1.
dist = @(x) 0;
for i = 1:length(surfFunc)
    dist = @(x) dist(x) + surfFunc(i)*x.^d;
    d = d - 1;
end
surfAileron = integral(dist, innerY, outerY);

CLSlope = resultInt/surfAileron/wantedAoA;

% still need to find Coeff -> lift from aileron sections (LLWING)
rollSlope = 2*CLSlope*eff/(wingS*b)*intResult;
rollCoeff = rollSlope*meanAngle;
totalRM = q*wingS*rollCoeff*b;
rotAngAcc = (totalRM/Ixx);


%% Pss and P calculus
% PROBLEM: STEADY STATE IS REALLY FAR. DRAG CAN'T EQUAL THE MOMENTUM.
% ROLLING IS AN ACCELERATING MOVEMENT.
CDr = 1.2;
fracHTail = 0.22;
fracVTail = 0.2;
S = wingS*(1 + fracHTail + fracVTail);
yD = 6; % conservative value.

% Pss reached when rolling moment equals drag rolling moment.
Pss = sqrt(2*totalRM/(rho*S*CDr*yD^3));     % steady state angular velocity

% integral 12.32 result for 0 - Pss 
SSAngle = (Ixx/(rho*S*CDr*yD^3))*log(Pss^2);

% simplification made:  tss not accounted for. 
PssDer = Pss^2/(2*SSAngle);
tss = sqrt(2*SSAngle/PssDer);

rollAngle = @(t) SSAngle + Pss*(t-tss);