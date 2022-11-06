clc;
clear; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Aileron rolling moment calculations

g = 9.81;
b747MTOW = 351500*g; % https://www.boeing.com/resources/boeingdotcom/commercial/airports/acaps/747_123sp.pdf
b747xx = 24675886;
MTOW = 39850*9.81;
Ixx = b747xx*(MTOW/b747MTOW);

% Wing characteristics
b = 11.85*2;
wingS = 65.258;
Cr = 4.08;
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
yo = innerEdge*b/2;
chordFInt = @(y) y^2/2 + 2/3*(lambda - 1)/b*y^3;
intResult = chordFInt(yo) - chordFInt(yi);

%% Aileron lift sections



% still need to find Coeff -> lift from aileron sections (LLWING)
rollSlope = 2*Coeff*eff/(wingS*b)*intResult;

rollCoeff = rollSlope*meanAngle;

totalRM = q*wingS*rollCoeff*b;

rotSpeed = (totalRM/Ixx)/(2*pi);


