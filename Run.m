close all
clear all
addpath('Functions')

% Initialize necessary parameters
% Default values and descriptions found in 'Functions/Initialize.m'
Parm = Initialize();

% Defines target output angle
target_angle = 50;

% Device parameters
Parm.Input.Wavelength = 1050;
Parm.Input.Polarization = 'TM';
Parm.Optimization.Target = [1];
Parm.Geometry.Thickness = 325; % Device layer thickness

% Compute necessary period corresponding to target angle
period = [Parm.Input.Wavelength*Parm.Optimization.Target/(sind(target_angle)-sind(Parm.Input.Theta)),0.5*Parm.Input.Wavelength];
Parm.Geometry.Period = period;

% Define # of Fourier orders
Parm.Simulation.Fourier = [12 12];

% Plot efficiency history
Parm.Display.PlotEfficiency = 1;

% Run one-step simulation
out = OneSimulation(Parm)
