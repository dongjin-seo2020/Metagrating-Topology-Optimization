% Single wavelength deflector optimization
% See Initialize() for definitions of parameters in Parm
function Out = OneSimulation(Parm)

%Extract common values for easier use
Wavelength = Parm.Input.Wavelength;
Period = Parm.Geometry.Period;
nBot = Parm.Geometry.Substrate;
nTop = Parm.Geometry.Top;
nDevice = Parm.Geometry.Device;

%Compute incident k-vector
kParallelForward = sind(Parm.Input.Theta);

% Compute total Fourier orders
NFourier = ceil(Parm.Simulation.Fourier.*Period/Wavelength);

% Define polarization values
if strcmp(Parm.Input.Polarization,'TE') 
    Polarizations = 1;
elseif strcmp(Parm.Input.Polarization,'TM')
    Polarizations = -1;
elseif strcmp(Parm.Input.Polarization,'Both')
    Polarizations = [1, -1];
else
    error('Invalid polarization');
end
NumPol = length(Polarizations); 

% Define grid for the device
[xGrid, yGrid, GridScale] = DefineGrid(Parm.Simulation.Grid, Period, Wavelength);
Nx = length(xGrid); %Number of x grid points
Ny = length(yGrid); %Number of y grid points

% If no starting point is given, generate a random starting point
if isempty(Parm.Optimization.Start)
    DevicePattern = RandomStart(Nx,Ny,Period,...
        Parm.Optimization.RandomStart,Parm.Geometry.SymmetryX,Parm.Geometry.SymmetryY);
else % Else regrid and use given starting point
    DeviceIn = Parm.Optimization.Start;
    DevicePattern = FineGrid(DeviceIn,Period,[Nx Ny]./size(DeviceIn),0,1); 
end
StartPattern = DevicePattern;

% Define full device stack
DeviceProfile = {[0,Parm.Geometry.Thickness,0],[1,3,2]}; % See Reticolo documentaion for definitions

% Generate binarization parameter B
BVector = GenerateBVector(Parm.Optimization.Binarize);

% Generate thresholding parameters for robustness
[ThresholdVectors, NRobustness] = GenerateThreshVectors(Parm);

% Compute blur radii in grid units
BlurGridLarge = Parm.Optimization.Filter.BlurRadiusLarge/GridScale;
BlurGrid = Parm.Optimization.Filter.BlurRadius/GridScale;
Figs = [];
% Initializing plot for geometries
if Parm.Display.PlotGeometry 
    Figs.FigGeo = figure;
end

% Initializing plot for geometries
if Parm.Display.PlotEfficiency
    Figs.FigEff = figure;
end

% Store efficiency at each iteration
AbsoluteEfficiency = zeros(NRobustness,NumPol);
RelativeEfficiency = zeros(NRobustness,NumPol);

iterStart = 1;
% Load checkpoint file if exists
if Parm.Checkpoint.Enable
    CheckpointFile = Parm.Checkpoint.File;
    if exist(CheckpointFile, 'file')
        load(CheckpointFile);
        iterStart = iter;
    end
end

%Initialize Reticolo
retio([],inf*1i);

%Main optimization loop
for iter = iterStart:MaxIterations
    
    % Save checkpoint file if appropriate
    if Parm.Checkpoint.Enable && (mod(iter, Parm.Checkpoint.Frequency) == 0)
        save(CheckpointFile);
    end
    tic;
    
    % First filter to enforce binarization
    FilteredPattern = DensityFilter2D(DevicePattern,BlurGridLarge);
    BinaryPattern = ThreshFilter(FilteredPattern,BVector(iter),0.5);
  
    GradientsAll = cell([NRobustness, 1]);
    DispPattern = cell([NRobustness, 1]);
    StoredAbsoluteEfficiencies = zeros(NRobustness, NumPol)
    
    
    % Begin robustness loop
    % Can be changed to parfor as necessary
    for robustIter = 1:NRobustness
        % Second filter to model physical edge deviations
        FilteredPattern2 = GaussFilter2D(BinaryPattern,BlurGrid);
        FinalPattern = ThreshFilter(FilteredPattern2,BVector(iter),ThresholdVectors(robustIter, iter));
        DispPattern{robustIter} = FinalPattern;
        % Define textures for each layer
        LayerTextures = cell(1,3);
        LayerTextures{1} = {nTop};
        LayerTextures{2} = {nBot};
        nPattern = FinalPattern*(nDevice - nTop) + nTop;
        LayerTextures{3} = FractureGeom(nPattern,nTop,nDevice,xGrid,yGrid);

        % Initialize empty field matrix
        FieldProductWeighted = zeros(NumPol,Parm.Simulation.ZGrid,Nx,Ny);
        
        % Begin polarization loop
        % Can be changed to parfor as necessary
        for polIter = 1:NumPol  
            % Set simulation parameters in Reticolo
            ReticoloParm = SetReticoloParm(Parm, Polarizations, polIter);

            % res1 computes the scattering matrices of each layer
            LayerResults = res1(Wavelength,Period,LayerTextures,NFourier,kParallelForward,0,ReticoloParm);

            % res2 computes the scattering matrix of the full device stack
            DeviceResults = res2(LayerResults,DeviceProfile);

            if Parm.Optimization.Target(polIter) == 0 && Parm.Input.Theta == 0
                FieldConvention = -1; % For normal output, opposite sign convention
            else
                FieldConvention = 1;
            end

            if (Polarizations(polIter)==1) %For TE polarization
                % Extract simulation results
                TransmittedLight = DeviceResults.TEinc_bottom_transmitted;
                
                % Find appropriate target
                TargetIndex = find((TransmittedLight.order(:,1)==Parm.Optimization.Target(polIter))&(TransmittedLight.order(:,2)==0));
                
                % Store efficiencies
                AbsEff = TransmittedLight.efficiency_TE(TargetIndex);
                RelativeEfficiency(robustIter,polIter) = TransmittedLight.efficiency_TE(TargetIndex)/sum(sum(TransmittedLight.efficiency));
    
                % Compute input field for adjoint simulation
                AdjointIncidence = [0,FieldConvention*exp(1i*angle(conj(TransmittedLight.amplitude_TE(TargetIndex))))];
                normalization = sqrt(2/3); % Normalize field for polarizations

                % res3 computes internal fields for each layer
                [ForwardField,~,~] = res3(xGrid,yGrid,LayerResults,DeviceProfile,[0,1],ReticoloParm);
                
            elseif (Polarizations(polIter)==-1) %For TM polarization
                % Extract simulation results
                TransmittedLight = DeviceResults.TMinc_bottom_transmitted;
                
                % Find appropriate target and store efficiencies
                TargetIndex = find((TransmittedLight.order(:,1)==Parm.Optimization.Target(polIter))&(TransmittedLight.order(:,2)==0));
                AbsEff = TransmittedLight.efficiency_TM(TargetIndex);
                RelativeEfficiency(robustIter,polIter) = TransmittedLight.efficiency_TM(TargetIndex)/sum(sum(TransmittedLight.efficiency));
                
                % Compute input field for adjoint simulation
                AdjointIncidence = [-FieldConvention*exp(1i*angle(conj(TransmittedLight.amplitude_TM(TargetIndex)))),0];
                normalization = (3/2)*sqrt(2/3); % Normalize field for polarizations

                % res3 computes internal fields for each layer
                [ForwardField,~,~] = res3(xGrid,yGrid,LayerResults,DeviceProfile,[1,0],ReticoloParm);
            end

            % Sum gradient over all robustness variants
            Gradients=zeros(size(DevicePattern));
            for robustIter = 1:NRobustness
                Gradients = Gradients + GradientsAll{robustIter};
            end
            Gradient = Gradients;

            % Output the calculated values
            disp('Adjoint Gradient:');
            disp(Gradient);
            disp('Absolute Efficiency:');
            disp(AbsoluteEfficiency);
            disp('Relative Efficiency:');
            disp(RelativeEfficiency);

             % Save outputs
            Out.BasePattern = DevicePattern;
            Out.BinaryPattern = FilteredPattern;
            Out.FinalPattern = FinalPattern;
            Out.StartPattern = StartPattern;
            Out.AbsoluteEfficiency = AbsoluteEfficiency;
            Out.RelativeEfficiency = RelativeEfficiency;
end
