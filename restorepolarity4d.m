function [vOutput, vIndex] = restorepolarity4d(v, mode)
% RESTOREPOLARITY4D Restore polarity of 4D IR MRI image by inverting all
% positive values before the signal minimum.
%
%    [vOutput, vIndex] = restorepolarity4d(v, mode)
%
%    INPUT:
%
%    v       - 4D matrix with magnitude data, diemnsions organized as
%              x, y, z, TI.
%    mode    - Char representing wheter the mnimum should be inverted.
%
%              'includeMinimum' - Invert the minimum value.
%              'excludeMinimum' - Do not invert the minimum value.
%
%    OUTPUT:
%
%    vOutput - 4D matrix with magnitude data with restored polarity.
%    vIndex  - 3D index matrix representning the position of minimum in the
%              original data.

% Reshape input 4D matrix into 2D matrix (x * y * z, TI)
vLong = reshape(v, [size(v, 1) * size(v, 2) * size(v, 3), size(v, 4)]);
nVoxel = size(vLong, 1);

% Restore polarity
vLongOutput = zeros(size(vLong));
vLongIndex = zeros(size(vLong));
parfor iVoxel = 1 : nVoxel

    [restoredSignal, signalFlipIndex] = restorepolarity(vLong(iVoxel, :), mode);
    vLongOutput(iVoxel, :) = restoredSignal;
    vLongIndex(iVoxel, :) = signalFlipIndex;

end
vOutput = reshape(vLongOutput, [size(v, 1)  size(v, 2)  size(v, 3) size(v, 4)]);
vIndex = reshape(vLongIndex, [size(v, 1)  size(v, 2)  size(v, 3) size(v, 4)]);

end

%
% --- Functions ---
%

function [restoredSignal, signalFlipIndex] = restorepolarity(signal, mode)
% RESTOREPOLARITY Restore polarity of a single voxel by inverting all
% positive values before the signal minimum.
%
%    [restoredSignal, signalFlipIndex] = restorepolarity(signal, mode)
%
%    INPUT:
%
%    signal           - vector with single voxel signal.
%    mode             - string representing wheter the mnimum should be
%                       inverted.
%
%                       'includeMinimum' - Invert the minimum value.
%                       'excludeMinimum' - Do not invert the minimum value.
%
%    OUTPUT:
%
%    restoredSignal   - vector with signal with restored polarity.
%    signalFlipIndex  - index representning the position of minimum in the
%                       original signal.


restoredSignal = signal;

switch mode
    case 'includeMinimum'
        [~, signalFlipIndex] = min(signal);
        restoredSignal(1:signalFlipIndex) = restoredSignal(1:signalFlipIndex) * -1;
    case 'excludeMinimum'
        [~, signalFlipIndex] = min(signal);
        signalFlipIndex = signalFlipIndex - 1;
        restoredSignal(1:signalFlipIndex) = restoredSignal(1:signalFlipIndex) * -1;
    otherwise
        warning('Unknow mode specified!')
end

end