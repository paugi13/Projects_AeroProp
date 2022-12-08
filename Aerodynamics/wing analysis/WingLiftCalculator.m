%% Wing lift calculus
% Main goal: to calculate stall speed
clc; clear; close all;
format long;
addpath(genpath(fileparts(mfilename('fullpath'))));

opt = input('Add flap to wing? Y/N (1/0): ');
wantedAoA = input('Wing incidence angle?: ');

numSt = buildStringAD(wantedAoA);
if opt == 0
    direct = join(['wing analysis/workspaces/wingLiftdist', numSt]);
    load(direct);
else
    direct = join(['wing analysis/workspaces/wingLiftdistFlap', numSt]);
    load(direct);
end

rho = 1.225;
MTOW = 39850*9.81;
S = 65.258;

stallSpeed = sqrt(MTOW/(0.5*rho*S*wingCL));     % [m/s]
stallSpeed = stallSpeed*3.6;    % [km/h]

vLanding = 196.2128/3.6;
q = 0.5*vLanding^2*rho;
CLMAx = MTOW/(q*S);

% 2.5784
CLMAXFlap = 3.2924;
CLMAXSimple = 1.907;

ClMAXFlap = 3.6;
ClMAXSimple = 1.8;

deltaCl = ClMAXFlap - ClMAXSimple;
deltaCL = CLMAXFlap - CLMAXSimple;

Sf = deltaCL/(0.92*cosd(22.5)*deltaCl)*S;
