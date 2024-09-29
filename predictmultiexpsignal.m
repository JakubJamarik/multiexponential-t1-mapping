function [predictedSignalMap] = predictmultiexpsignal(t1Map, wMap, ti)
% PREDICTMUTIEXPSIGNAL Predicts a  multi-exponential signal based on T1 and
% W parameter maps. The signal is computed according to the following
% equation:
%
% S(TI) = w1 * (1 - 2 * exp( - TI / T1_1)) + w2 * (1 - 2 * exp( - TI /
% T1_2)) + ...
%
% where w represents the realtive weight of each component.
%
% [predictedSignalMap] = predictmultiexpsignal(t1Map, wMap, ti)
%
%    INPUT:
%
%    t1Map              - 4D array [x, y, z, numComponents], T1 times (ms)
%                         per voxel component.
%    wMap               - 4D array [x, y, z, numComponents], weights for
%                         each voxel component.
%    ti                 - Vector with inversion recovery times in ms.
%
%    OUTPUT:
%
%    predictedSignalMap -  4D array [x, y, z, numTimePoints] with the
%                          predicted signal.

%
% Main
%

% Get the size of the input arrays
[xDim, yDim, zDim, ~] = size(t1Map);
nTimePoint = length(ti);

% Pre-allocate the output array for the simulated signal
predictedSignalMap = zeros(xDim, yDim, zDim, nTimePoint);

% Loop over time points to compute the signal for each time point
for iTimePoint = 1:nTimePoint
    % Extract the current inversion time
    iTi = ti(iTimePoint);
    
    % Compute the signal for all components at the current time point
    % This calculation is performed element-wise across all voxels and components
    signalPerComponentMap = wMap .* (1 - 2 * exp(-iTi ./ t1Map));
    
    % Sum over the components dimension (4th dimension) to get the total signal
    predictedSignalMap(:, :, :, iTimePoint) = sum(signalPerComponentMap, 4);
    
end