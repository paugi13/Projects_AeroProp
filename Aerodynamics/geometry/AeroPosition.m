% Integration script for aircraft geometry 
% Script to get aircraft's geometry data for the integration part. 

clc; clear; close all;
format long;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Load CP coordinates. 
load('wing analysis/workspaces/WingCPCoords');
load('wing analysis/workspaces/HTailCPCoords');
load('wing analysis/workspaces/VTailCPCoords');

rootChord = 4.08;
% CP distances to wing's CP.
lh = 15.01;
lv = 14.2;

FrontCAWing = 13.9085;
%% Calculations for distances with respect to the wing's leading edge.
WingHorTE = lh - HCPCoords(2) + CPCoords(2) - rootChord;
WingVerTE = lv - VCPCoords(2) + CPCoords(2) - rootChord;

WingHorLE = lh - HCPCoords(2) + CPCoords(2);
WingVerLE = lv - VCPCoords(2) + CPCoords(2);

HTailEndCoord = FrontCAWing + lh + (HrootChord - HCPCoords(2));
VTailEndCoord = FrontCAWing + lv + (VrootChord - VCPCoords(2));