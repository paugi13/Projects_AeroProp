function [ResMoment, TailS, CMalpha] = PitchingEffectiveness(XcgMTOW, WingTailD,Xwing, XtailCG, ... 
    flightReg, tailIncidence, tailClSlope, VolumeCoeff, FuselageAoA, TailEff)
% Function to calculate resultant pitching moment taking geometry and tail
% properties as inputs.
WingS = 65.258;     % wing's surface
MAC = 2.986225279129153;     % mean aerodynamic chord 
AR = 8.6;
WingIncidence = 5.25;
TailS = VolumeCoeff*MAC*WingS/WingTailD;

% Load wing data for selected regime
if flightReg == 1   %takeoff (flap)
    ClSlope = [0.0746 1.66];    % alpha functions
    Cmac = -0.179488921644864;
    TOSpeed = 270.392/3.6;
    rho = 0.974;
    q = 0.5*rho*TOSpeed^2;
else    %cruise (no flap)
    ClSlope = [0.0795 0.0816];      % alpha functions
    Cmac = -0.059817260761975;
    cruiseSpeed = 857/3.6;
    rho = 0.363918;
    q = 0.5*rho*cruiseSpeed^2;
end

% Wing contributions
wingCL = (ClSlope(1)*(FuselageAoA+WingIncidence) + ClSlope(2));
wingLift = q*WingS*wingCL;  %*(XcgMTOW-Xwing);
wingMoment = q*WingS*Cmac*MAC;

% Wing downwash
eps0 = 2*wingCL/(pi*AR);
epsAlphaSlope = 2*ClSlope(1)/(pi*AR);
eps = eps0 + epsAlphaSlope*FuselageAoA;

tailCL = tailClSlope(1)*(tailIncidence + FuselageAoA - eps) + tailClSlope(2);
tailLift = q*TailS*tailCL;
M1 = wingLift*(XcgMTOW-Xwing);
M2 = wingMoment;
M3 = tailLift*XtailCG;

% Pitching moment slope: < 0 means stable conf.
CMalpha = ClSlope(1)*(XcgMTOW-Xwing) + tailClSlope(1)*(TailS/WingS)*XtailCG*(1-...
    epsAlphaSlope + TailEff*DE_Flap);

% Resultant moment
ResMoment = M1 + M2 + M3;
end

