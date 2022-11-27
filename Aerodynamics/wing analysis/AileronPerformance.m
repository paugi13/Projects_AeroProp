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
rootThick = 0.14;
tAverage = Cr*rootThick*0.5;
Ixx = MFuselage*R^2 + 1/12*MWing*(b^2 + tAverage^2);

% Wing characteristics
wingS = 65.258;
lambda = 0.35;
chordF = @(y) Cr*(1 + 2*(lambda - 1)/b*y); 
Ct = chordF(b/2);
upAngle = 20; 
downAngle = 20; 
meanAngle = (abs(upAngle) + abs(downAngle))/2;

% Aileron characteristics
xh = 0.7;
eff = 0.6566; %0.6566;
innerEdge = 0.7;
outerEdge = 0.95;

% Flight profile data
turnSpeed = 350/3.6;
rho = 0.974;
q = 0.5*rho*turnSpeed^2;

%% Rolling moment integration 
yi = innerEdge*b/2;
yo = outerEdge*b/2;
chordFInt = @(y) Cr*(y^2/2 + 2/3*(lambda - 1)/b*y^3);
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
rotAngAcc = (totalRM/Ixx)*180/pi;
requiredTime30 = sqrt(30*2/rotAngAcc);

disp([num2str(requiredTime30), ' s']);
disp([num2str(rotAngAcc), ' degress per squared second']);
