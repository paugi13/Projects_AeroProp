% Integration script for aircraft geometry 
% Script to get aircraft's geometry data for the integration part. 

clc; clear; close all;
format long;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Load CP coordinates. 
load('wing analysis/workspaces/WingCPCoords');
load('wing analysis/workspaces/HTailCPCoords');
load('wing analysis/workspaces/VTailCPCoords');

% CP distances to wing's CP.
lh = 12.991;
lv = 12.2;

%% Calculations for distances with respect to the wing's leading edge.
WingHor = lh - HCPCoords(2) + CPCoords(2);
WingVer = lv - VCPCoords(2) + CPCoords(2);