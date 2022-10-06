
% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodinàmica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

clc; clear all; close all;
format long;

% -------------------------------------------------------------------------
% INPUT DATA
% -------------------------------------------------------------------------

% Wing planform (assumes planar wing)

AR = 21.3 ;   % aspect ratio
TR = 1/5.55 ;   % taper ratio
DE25 = 17 ; % sweep angle at c/4 (deg)

ETIP = -7.1; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip) -->
% AIRFOIL DATA

A0p = [ -1 -1 ]; % root and tip section zero-lift angles (deg)
CM0p = [ 0.012 0.012 ]; % root and tip section free moments
CDP = [ 0.0066 -0.0061 0.0063 ;   % root section CD0, k1 and k2  (airfoil CD curve)
         0.01 -0.0076 0.0096 ] ;  % tip section CD0, k1 and k2

% Flap/aileron (symmetrical deflection)

YF_pos = [ 0.0 0.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2 ;  % flap_chord/chord ratio
DE_flap = 0.0; % flap deflection (deg, positive:down)
FlapCorr = 1.0 ; % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)

N = 100 ; % number of panels along the span

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ] ; % angles of attack for analysis (deg) 

% -------------------------------------------------------------------------
% LIFTING LINE SOLUTION
% -------------------------------------------------------------------------

% Wing discretization (lenghts are dimensionless with the wing span)

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

% Assembly of the influence coefficients matrix (needed once)

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

% Solve circulations for different angles of attack

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

% Loads calculation using plane Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

% -------------------------------------------------------------------------
% POSTPROCESS ( your work starts here...)
% -------------------------------------------------------------------------

