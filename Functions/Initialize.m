% Initialize default values for necessary optimization parameters
function [Parm] = Initialize()
addpath('reticolo_allege')
Parm.Input.Wavelength = 1000; % Wavelength of incident light
Parm.Input.Polarization = 'TM'; % Polarization: can be 'TM', 'TE' or 'Both'
Parm.Input.Theta = 0; % Incident angle (angle given in air)

Parm.Geometry.Period = [1500, 500]; % Period of the device
Parm.Geometry.Thickness = 350; % Device layer thickness
Parm.Geometry.Substrate = 1.45; % Refractive index of substrate
Parm.Geometry.Device = 3.45; % Refractive index of device
Parm.Geometry.Top = 1; % Refractive index of upper layer

% Symmetry enables faster computation
Parm.Geometry.SymmetryY = 1; % Enforces symmetry in the Y direction
Parm.Geometry.SymmetryX = 0; % Enforces symmetry in the X direction

Parm.Simulation.Grid = 90; % Grid points per wavelength 
Parm.Simulation.Fourier = [14 14]; % Number of Fourier orders per wavelength in x and y respectively
Parm.Simulation.ZGrid = 50; % Number of slices in device layer at which fields are computed

Parm.Optimization.Target = 1; % Target diffraction order in the x-direction for each polarization i.e. [1 -1]

Parm.Optimization.Iterations = 350; % Max iterations 
Parm.Optimization.Gradient.StepSize = 0.1; % Initial gradient step size
Parm.Optimization.Gradient.StepDecline = 0.992; % Multiplying factor that decreases step size each iteration

Parm.Optimization.Filter.BlurRadiusSmall = 5; % Blur radius for small blurring
Parm.Optimization.Filter.BlurRadius = 20; % Blurring radius used for second blurring
Parm.Optimization.Filter.BlurRadiusLarge = 60; % Blurring radius for initial blurring
Parm.Optimization.Filter.BlurSmallIter = 5; % Frequency of small blurring
Parm.Optimization.Filter.BlurSmallIterStop = 25; % Steps at the end of optimization which have no blurring
Parm.Optimization.Filter.BlurLargeIter = 60; % Frequency of large blurring
Parm.Optimization.Filter.BlurLargeIterStop = 50; % Steps at the end of optimization which have no blurring

% Parameters for defining the binarization rate B. Smaller B is a
% less-sharp function, a large B will create a step function that cuts off 
% at the different values of N.
Parm.Optimization.Binarize.Min = 1; % Minimum value of B that is used after it is initialized. 
Parm.Optimization.Binarize.Max = 50; % Maximum value of B
Parm.Optimization.Binarize.IterationStart = 1; % Iteration at which B is allowed to be nonzero.
Parm.Optimization.Binarize.IterationHold = 20; % Number of iterations during which B is kept constant before increasing (after it is initialized)

% Start and End should be the same length. Each element represents a
% separate thread of the optimization at a different flat edge deviation in nm. 
Parm.Optimization.Robustness.StartDeviation = [-5 0 5]; % Starting edge deviation values
Parm.Optimization.Robustness.EndDeviation = [-5 0 5]; % Ending edge deviation values
Parm.Optimization.Robustness.Ramp = 30; % Iterations over which the thresholding parameter changes
Parm.Optimization.Robustness.Weights = [.5 1 .5]; % Gradient weight for each robustness value

Parm.Optimization.Start = []; % Input pattern matrix as starting point. Leave empty to generate random starting point

% Parameters for random starting point. See RandomStart for details. 
% Ignored if starting point is already provided
Parm.Optimization.RandomStart.Pitch = 100; % Grid size on which to create random Gaussians
Parm.Optimization.RandomStart.Average = 0.6; % Average refractive index
Parm.Optimization.RandomStart.Sigma = 0.5; % Standard deviation of the refractive index

Parm.Display.ShowText = 1; % Display progress as console output
Parm.Display.PlotGeometry = 1; % Display current device geometry as image
Parm.Display.PlotEfficiency = 0; % Plot device efficiency vs iteration

% Checkpointing saves the optimization state to a temporary file 
% and allows restarting from the previous checkpoint if the optimization is
% interrupted
Parm.Checkpoint.Enable = false; % Enable checkpointing behavior
Parm.Checkpoint.File = ''; % Checkpoint file name
Parm.Checkpoint.Frequency = 10; % Number of iterations between checkpoints

% Disable parfor warning messages
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
end