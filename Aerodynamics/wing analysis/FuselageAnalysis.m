% Code for fuselage aerodynamic performance analysis.
clc;
clear;
close all;

WantedAoA = input('Wanted cruise AoA? ');


%% Flight data
ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10.] ;

for i = 1:length(ALPHA)
    if ALPHA(i) == WantedAoA
        pos = i;
        a = 1;
        break
    end
end

if a~=1
    error('Angle of attack not in the alpha vector');
end

%%  Drag analysis
% Fuselage dimensions.
r = 1.4293;
b = 2*r;
lh = 31.26320;

Re=4e6;
Cf=0.074/(Re^(1/5));

Ss = 219.62;    %   Total wetted area (ellipsoid area)
S_bal = pi*r.*sqrt((1+(tand(ALPHA)).^2).*(1./(1/r^2+(tand(ALPHA)).^2)));
% Frontal area as a function of AoA.

Cd0 = Cf*(1+60/((lh/b)^3)+0.0025*(lh/b))*(Ss./S_bal);

%% Pitching moment coefficient
Vf = 133.76;
v = 800;
rho = 1.225;
q = 0.5*rho*v^2;
k = 1;

Cmf = (pi*k*q*2*0.09*ALPHA*pi/180*4/3)/(q*Vf);

%% Drag and pitching moment computation. 
MfCruise = Cmf(pos)*q*Vf;
D = Cd0(pos)*q*(Vf^(2/3));
