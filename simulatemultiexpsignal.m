function simulatedSignalMap = simulatemultiexpsignal(gtT1Map, wMap, ti, snr)
% SIMULATEMULTIEXPSIGNAL Returns a simulated multi-exponential signal with
% added noise. The signal is computed according to the following equation:
%
% S(TI) = w1 * (1 - 2 * exp( - TI / T1_1)) + w2 * (1 - 2 * exp( - TI /
% T1_2)) + ...
%
% where w represents the realtive weight of each component.
%
%    [simulatedSignalMap] = simulatemultiexpsignal(gtT1Map, wMap, ti, snr)
%
%    INPUT:
%
%    gtT1Map            - 4D array [x, y, z, numComponents], T1 times (ms)
%                         per voxel component.
%    wMap               - 4D array [x, y, z, numComponents], weights for
%                         each voxel component.
%    ti                 - Vector with inversion recovery times in ms.
%    snr                - Signal-to-noise ratio in dB.
%
%    OUTPUT:
%
%    simulatedSignalMap -  4D array [x, y, z, numTimePoints] with the
%                          simulated signal.

%
% Main
%

% Get the size of the input arrays
[xDim, yDim, zDim, ~] = size(gtT1Map);
nTimePoint = length(ti);

% Pre-allocate the output array for the simulated signal
simulatedSignalMap = zeros(xDim, yDim, zDim, nTimePoint);

% Loop over time points to compute the signal for each time point
for iTimePoint = 1:nTimePoint
    % Extract the current inversion time
    iTi = ti(iTimePoint);
    
    % Compute the signal for all components at the current time point
    % This calculation is performed element-wise across all voxels and components
    signalPerComponentMap = wMap .* (1 - 2 * exp(-iTi ./ gtT1Map));
    
    % Sum over the components dimension (4th dimension) to get the total signal
    simulatedSignalMap(:, :, :, iTimePoint) = sum(signalPerComponentMap, 4);
end

% Add Gaussian noise based on SNR
% Convert SNR from dB to linear scale
snrLinear = 10^(snr / 20);

% Calculate the standard deviation of the noise from SNR
signalPower = mean(abs(simulatedSignalMap(:)).^2);  
noiseStd = sqrt(signalPower) / snrLinear;

% Add Gaussian noise with mean 0 and standard deviation `noiseStd`
noise = noiseStd * randn(size(simulatedSignalMap));
simulatedSignalMap = simulatedSignalMap + noise;

end
