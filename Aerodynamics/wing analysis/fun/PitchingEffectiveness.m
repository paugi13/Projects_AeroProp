function [M1, M2, M3, ResMoment, TailS, CMalpha, eps] = PitchingEffectiveness(XcgMTOW, ...
    WingTailD,Xwing, XtailCG, flightReg, tailIncidence, tailClSlope, ...
    VolumeCoeff, FuselageAoA, TailEff, DE_Flap)
% Function to calculate resultant pitching moment taking geometry and tail
% properties as inputs.
WingS = 65.258;     % wing's surface
MAC = 2.986225279129153;     % mean aerodynamic chord 
AR = 8.6;
WingIncidence = 5.25;
TailS = VolumeCoeff*MAC*WingS/WingTailD;
mu = 0.97;

% Load wing data for selected regime
if flightReg == 1   %takeoff (flap)
    ClSlope = [0.0746 1.66];    % alpha functions
    Cmac = -0.179488921644864;
    TOSpeed = 270.392/3.6;
    rho = 0.974;
    q = 0.5*rho*TOSpeed^2;
else    %cruise (no flap)
    ClSlope = [0.0795 0.0816];      % alpha functions
    Cmac = -0.062685;                  %0.059817260761975; ;
    cruiseSpeed = 857/3.6;
    rho = 0.363918;
    q = 0.5*rho*cruiseSpeed^2;
end

% Wing contributions
wingCL = (ClSlope(1)*(FuselageAoA+WingIncidence) + ClSlope(2));
wingLift = q*WingS*wingCL;  %*(XcgMTOW-Xwing);
wingMoment = q*WingS*Cmac*MAC;

% Wing downwash
eps0 = 2*ClSlope(2)/(pi*AR);
epsAlphaSlope = 2*ClSlope(1)*180/pi/(pi*AR);
eps = (eps0 + epsAlphaSlope*((FuselageAoA + WingIncidence)*pi/180))*180/pi;

tailCL = tailClSlope(1)*(tailIncidence + FuselageAoA - eps + ...
    TailEff*DE_Flap) + tailClSlope(2);
tailLift = q*TailS*tailCL*mu;
M1 = wingLift*(XcgMTOW-Xwing);
M2 = wingMoment;
M3 = tailLift*XtailCG;

% Pitching moment slope: < 0 means stable conf.
% CMalpha = ClSlope(1)*(XcgMTOW-Xwing) + tailClSlope(1)*(TailS/WingS)*XtailCG*(1-...
%     epsAlphaSlope);
CMalpha = ClSlope(1)*(XcgMTOW-Xwing) - tailClSlope(1)*(1-epsAlphaSlope)*...
    VolumeCoeff*0.95;

% T1 = ClSlope(1)*(XcgMTOW-Xwing);
% T2 = tailClSlope(1)*(TailS/WingS)*XtailCG*(1-...
%     epsAlphaSlope)/MAC;
% Resultant moment
ResMoment = M1 + M2 + M3;
end

