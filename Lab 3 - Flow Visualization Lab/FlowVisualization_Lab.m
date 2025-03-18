%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reynaldo Villarreal Zambrano
% February 5, 2025
% TFES Lab 3 - Flow Visualization Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all

%% Declare constants

% kinematic viscosity
mu = 1.0 * 10^-6; % m^2 / s (@ 8 degrees C, water)

% declare function to find flow speed given frequency (x)
flow_speed = @(x) 0.175 + 1.952 * x; % cm/s

%% calculate flow speed for circular cylinder
cylinder_dim = 13; % diameter (mm)
cylinder_dim = cylinder_dim / 1000; % mm -> m
cylinder_low_speed = flow_speed(4); % 4 Hz -> cm/s
cylinder_low_speed = cylinder_low_speed / 100; % cm/s -> m/s
cylinder_high_speed = flow_speed(8); % 8 Hz -> cm/s
cylinder_high_speed = cylinder_high_speed / 100; % cm/s -> m/s

% find reynold's number for cylinder flows
cylinder_low_reynolds = (cylinder_low_speed * cylinder_dim) / mu;
cylinder_high_reynolds = (cylinder_high_speed * cylinder_dim) / mu;

% print results
fprintf("Reynold's number for low speed around a cylinder: %.4f\n", cylinder_low_reynolds);
fprintf("Reynold's number for high speed around a cylinder: %.4f\n", cylinder_high_reynolds);


%% calculate flow speed for long flat plate
plate_dim = 455; % length (mm)
plate_dim = plate_dim / 1000; % mm -> m
plate_low_speed = flow_speed(4); % 4 Hz -> cm/s
plate_low_speed = plate_low_speed / 100; % cm/s -> m/s
plate_med_speed = flow_speed(8); % 8 Hz -> cm/s
plate_med_speed = plate_med_speed / 100; % cm/s -> m/s
plate_high_speed = flow_speed(12); % 12 Hz -> cm/s
plate_high_speed = plate_high_speed / 100; % cm/s -> m/s

% calculate reynold's number for flat plate flow
plate_low_reynolds = (plate_low_speed * plate_dim) / mu; 
plate_medium_reynolds = (plate_med_speed * plate_dim) / mu;
plate_high_reynolds = (plate_high_speed * plate_dim) / mu;

fprintf("Reynold's number for medium speed over a flat plate: %.4f\n", plate_medium_reynolds);
fprintf("Reynold's number for low speed over a flat plate: %.4f\n", plate_low_reynolds);
fprintf("Reynold's number for high speed over a flat plate: %.4f\n", plate_high_reynolds);

%% calculate flow speed for airfoil
airfoil_dim = 50; % width (mm)
airfoil_dim = airfoil_dim / 1000; % mm -> m
airfoil_speed = flow_speed(8); % 8 Hz -> cm/s
airfoil_speed = airfoil_speed / 100; % cm/s -> m/s

% calculate reynolds number for airfoil
airfoil_reynolds = (airfoil_speed * airfoil_dim) / mu;
fprintf("Reynold's number for same speed over airfoil: %.4f\n", airfoil_reynolds);

%% calculate flow speed for other shape (golf ball)
golfball_dim = 44; % diameter (mm)
golfball_dim = golfball_dim / 1000; % mm -> m
golfball_high_speed = flow_speed(12); % 12 Hz -> cm/s
golfball_high_speed = golfball_high_speed / 100; % cm/s -> m/s

% calculate reynolds number for golfball
golfball_high_reynolds = (golfball_high_speed * golfball_dim) / mu;
fprintf("Reynold's number for high speed over golfball: %.4f\n", golfball_high_reynolds);
