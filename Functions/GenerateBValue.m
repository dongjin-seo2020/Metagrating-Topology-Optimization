% Defines the values of B, the binarization rate, at each iteration
% Plot BVector to see the progression of binarization
function BVector = GenerateBValue(BinParm)
    % Extract parameters
    BMin = BinParm.Min;
    
    % Set the binarization rate to a constant value
    BValue = BMin;
end