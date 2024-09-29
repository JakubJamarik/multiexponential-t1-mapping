function [t1Map, nComponenMap] = calculatet1mapilt(v, ti, MIN_PEAK_HEIGHT, varargin)
% CALCULATET1MAPILT Returns the T1 map and the map of number of components
% of the input 4D IR MRI image, estiamted by the Discrete Inverse Laplace
% Transform. Only values of T1 with M0 higher than a treshold are retained.
%
%    [t1Map, nComponentMap] = calculatet1mapilt(v, ti, MIN_PEAK_HEIGHT, t1t, model)
%
%    INPUT:
%
% 	 v               - 4D (x, y, z, TI) matrix of IR MRI image.
% 	 ti              - Vector of TI times in ms.
%    MIN_PEAK_HEIGHT - Constat indicating the relative weight threshold for
%                      component retention.
%	 t1t             - Vector of T1 basis values. Default is
%                      linspace(50, 5000, 100)
%    model           - Char indicating which model to use. Default is
%                      'absoluteValueOfSum'.
%
%                      'sumOfAbsoluteValues' - The signal is modeled as a
%                                              sum of absolute values of
%                                              individual exponential
%                                              functions.
%                      'absoluteValueOfSum'  - The signal is modeled as an
%                                              absolute value of a sum of
%                                              individual exponential
%                                              functions.
%
%    OUTPUT:
%
%	 t1map           - 4D (x, y, z, t1t) parametric map of estimated T1
%                      values.
%    nComponenMap    - 4D (x, y, z, n) parametric map indicating number
%                      of components per voxel.

%
% Defaults
%

t1t = linspace(50, 5000, 100);
model = 'absoluteValueOfSum';
Defaults = {t1t, model};
Defaults(1:length(varargin)) = varargin;
[t1t, model] = Defaults{:};

%
% Main
%

% Compute T1 map
[t1Map, ~, normM0Map, ~] = calculatet1mapilt4D(v, ti, t1t, model);

% Treshold T1 map based on normM0Map (relative weight)
[t1Map, ~, nComponenMap] = tresholdt1map(t1Map, normM0Map, MIN_PEAK_HEIGHT);

% Drop T1s with zero normM0 across all voxels
t1Map = dropzerocomponents(t1Map, nComponenMap);

end

function [t1Map, m0Map, normM0Map, rmseMap] = calculatet1mapilt4D(v, ti, varargin)
% CALCULATET1MAPILT4D Estimate raw T1 and M0 maps using the
% Discrete Inverse Laplace Transform. This function is a wrapper for a
% single voxel computation.
%
%    [t1Map, m0Map, normM0Map, rmseMap] = calculatet1mapilt4D(v, ti, t1t, model)
%
%    INPUT:
%
% 	 v         - 4D (x, y, z, TI) matrix of IR MRI image.
% 	 ti        - Vector of TI times in ms.
%	 t1t       - Vector of T1 basis values.
%    model     - Char indicating which model to use. Default is
%                'absoluteValueOfSum'.
%
%                      'sumOfAbsoluteValues' - The signal is modeled as a
%                                              sum of absolute values of
%                                              individual exponential
%                                              functions.
%                      'absoluteValueOfSum'  - The signal is modeled as an
%                                              absolute value of a sum of
%                                              individual exponential
%                                              functions.
%
%    OUTPUT:
%
%	 T1map     - 4D (x, y, z, t1t) parametric map of estimated T1 values.
%	 M0map     - 4D (x, y, z, t1t) parametric map of estimated M0 values.
%    normM0Map - 4D (x, y, z, t1t) parametric map of estimated M0 values
%                normalized across the 4th dimension.
%    rmseMap   - 4D (x, y, z, t1t) parametric map of RMSE values.

%
% Defaults
%

t1t = linspace(50, 5000, 100);
model = 'absoluteValueOfSum';
Defaults = {t1t, model};
Defaults(1:length(varargin)) = varargin;
[t1t, model] = Defaults{:};

%
% Input validation
%

% Set input vectors as column vectors
if isrow(ti)
    ti = ti';
end
if isrow(t1t)
    t1t = t1t';
end

%
% Main
%

% Compute size of input volume
[xSize, ySize, zSize, tSize] = size(v);
t1tSize = length(t1t);

% Reshape input volume into 2D (x * y * z, TI)
vLong = reshape(v, [xSize * ySize * zSize, tSize]);

% Compute T1 and M0 maps
t1MapLong = zeros(size(vLong, 1), length(t1t));
m0MapLong = zeros(size(vLong, 1), length(t1t));
rmseMapLong = zeros(size(vLong, 1), 1);
parfor k = 1:(size(vLong, 1))
    
    % Progress counter
    if k == 1
        disp('Fit started ')
    elseif mod(k, round((size(vLong, 1)) * 0.1)) == 0
        disp('. ')
    elseif k == (size(vLong, 1))
        disp('100% done.')
    end
    
    % Estiamte M0 per voxel
    m0MapLong(k, :) = diltnnls(vLong(k, :)', ti, t1t, model);
    
    % Estimate correspondning T1
    t1MapLong(k, :) = (m0MapLong(k, :) > 0.0001)' .* t1t;
    
    % Compute RMSE
    r1MapLong = 1 ./ t1MapLong(k, :);
    kernelMapLong = - ti * r1MapLong;
    kernelMapLong = ones(size(kernelMapLong)) - 2 * exp(kernelMapLong);
    kernelMapLong = abs(kernelMapLong);
    signalEstimateMapLong = m0MapLong(k, :) * kernelMapLong';
    rmseMapLong(k, :) = sqrt(mean((vLong(k, :) - signalEstimateMapLong) .^ 2));
    
end

% Transform data back to 4D
t1Map = reshape(t1MapLong, [xSize ySize zSize t1tSize]);
m0Map = reshape(m0MapLong, [xSize ySize zSize t1tSize]);
rmseMap = reshape(rmseMapLong, [xSize ySize zSize 1]);

% Compute normalized M0 map
normM0Map = normalizematrix(m0Map, 4);

end

function [p] = diltnnls(signal, ti, varargin)
% DILTNNLS Estimate the weights (M0) for each basis function using the
% Discrete Inverse Laplace Transform implemented by Non Negative Least
% Squares.
%
%    [p] = diltnnls(signal, ti, t1t, model)
%
%    INPUT:
%
%    signal - Vector of IR signal dependent on TI.
%    ti     - Vector of TI times in ms.
%	 t1t    - Vector of T1 basis values.
%    model  - Char indicating which model to use. Default is
%             'absoluteValueOfSum'.
%
%             'sumOfAbsoluteValues' - The signal is modeled as a
%                                     sum of absolute values of
%                                     individual exponential
%                                     functions.
%             'absoluteValueOfSum'  - The signal is modeled as an
%                                     absolute value of a sum of
%                                     individual exponential
%                                     functions.
%
%    OUTPUT:
%
%    p      - Vector of weights (M0) for each basis T1 value.

%
% Defaults
%

t1t = linspace(50, 5000, 100);
model = 'absoluteValueOfSum';

% Overwrite defautls if supplied in varargin
Defaults = {t1t, model};
Defaults(1:length(varargin)) = varargin;
[t1t, model] = Defaults{:};

%
% Input validation
%

% Set input vectors as column vectors
if isrow(signal)
    signal = signal';
end
if isrow(ti)
    ti = ti';
end
if isrow(t1t)
    t1t = t1t';
end

%
%  Main
%

% Compute basis matrix
switch model
    case 'sumOfAbsoluteValues'
        t1t = 1 ./ t1t;
        K = - ti * t1t';
        K = ones(size(K)) - 2 * exp(K);
        K = abs(K);
    case 'absoluteValueOfSum'
        t1t = 1 ./ t1t;
        K = - ti * t1t';
        K = exp(K);
        K = ones(size(K)) - 2 * K;
    otherwise
        error("Unknow model!!!");
end

% Find regularization parameter
lambda = lcurvegcv(K, signal);

% Regularize and solve using quadprog
H = K' * K;
H = H + lambda * ones(size(H));
f = - K' * signal;
[~, n] = size(K);
p = quadprog(H, f, -eye(n,n), zeros(n,1));

%
% --- Functions ---
%

    function [lambda] = lcurvegcv(C, d)
        % LCURVEGCV Find the regularization parameter lambda by optimizaing the
        % L-curve via General Crossvalidation. The function is based on the
        % original matlab codewriten by Thomas Prasloski (email:
        % tprasloski@gmail.com). For details on the L-curve and GCV methods see
        % [Hansen, P.C., 1992. Analysis of Discrete Ill-Posed Problems by Means of the L-Curve. SIAM Review, 34(4), 561-580](https://doi.org/10.1137/1034115).
        %
        %    [lambda] = lcurvegcv(C, d)
        %
        %    INPUT:
        %
        %    C      - Basis matrix.
        %    d      - Vector of IR signal dependent on TI.
        %
        %    OUTPUT:
        %
        %    lambda - Regularization parameter.
        
        %
        % Main
        %
        
        lambda = fminbnd(@(x) G(x,C,d), 0, 0.1, optimset('TolX', 0.001));
        
        function [g] = G(lambda,C,d)
            % G Helper function to be optimized by lcurvegcv.
            %
            % [g] = G(lambda,C,d)
            %
            %    INPUT:
            %
            %    lambda - Initial value of regularization parameter.
            %    C      - Basis matrix.
            %    d      - Vector of IR signal dependent on TI.
            %
            %    OUTPUT:
            %
            %    g      - Output value of regularization parameter.
            
            g=sqrt(sum((d-C*lsqlin([C;lambda*eye(size(C,2))],[d;zeros(size(C,2),1)])).^2))/((trace(eye(size(C,1))-(C*((C'*C+lambda*eye(size(C,2)))^-1)*C')))^2);
            
        end
        
    end

end

function [t1MapTresholded, normM0MapTresholded, nPeakMap] = tresholdt1map(t1Map, normM0Map, MIN_PEAK_HEIGHT)
% TRESHOLDT1MAP Function returns T1 maps tresholded based on component
% relative weight (normalized M0 map). First, a set of peaks for each voxel
% is found, based on the treshold. All other components within the given
% voxel are then set to zero.
%
%    [t1MapTresholded, normM0MapTresholded, nPeakMap] = tresholdt1map(t1Map, normM0Map, MIN_PEAK_HEIGHT)
%
%    INPUT:
%
%	 T1map               - 4D (x, y, z, t1t) parametric map of estimated T1
%                          values.
%    normM0Map           - 4D (x, y, z, t1t) parametric map of estimated M0
%                          values normalized across the 4th dimension.
%    MIN_PEAK_HEIGHT     - Scalar indicating the treshold for realtive
%                          weight (M0).
%
%    OUTPUT:
%
%    t1MapTresholded     - 4D (x, y, z, t1t) parametric map of estimated T1
%                          values, tresholded based on their relative
%                          weight (M0).
%    normM0MapTresholded - 4D (x, y, z, t1t) parametric map of estimated M0
%                          values normalized across the 4th dimension,
%                          tresholded based on their relative weight (M0).
%    nPeakMap            - 3D (x, y, z) parametric map of the number of
%                          components.

% Get size of the input array.
[xSize, ySize, zSize, tSize] = size(t1Map);

% Reshape array into long format (2D matrix).
t1MapLong = reshape(t1Map, [xSize * ySize * zSize, tSize]);
normM0MapLong = reshape(normM0Map, [xSize * ySize * zSize, tSize]);
nVoxel = size(t1MapLong, 1);

t1MapLongOutput = t1MapLong;
normM0MapLongOutput = normM0MapLong;
for iVoxel = 1 : nVoxel
    
    % Act as if relative weights in normalized M0 map are signal 'peaks' and find their location based on threshold
    [~, peakLocation] = findpeaks(normM0MapLong(iVoxel, :), 'SortStr', 'descend', 'MinPeakHeight', MIN_PEAK_HEIGHT);
    
    peakLocationLength = size(peakLocation, 2);
    switch peakLocationLength
        case 0
            % If no peak was found set all estimates to zero
            t1MapLongOutput(iVoxel, :) = 0;
            normM0MapLongOutput(iVoxel, :) = 0;
        otherwise
            % If atleast one peak was found set all remaining estimates to zero
            t1MapLongOutput(iVoxel, setdiff(1:end, peakLocation)) = 0;
            normM0MapLongOutput(iVoxel, setdiff(1:end, peakLocation)) = 0;
    end
    
end

% Reshape output into original size.
t1MapTresholded = reshape(t1MapLongOutput, [xSize ySize zSize tSize]);
normM0MapTresholded = reshape(normM0MapLongOutput, [xSize ySize zSize tSize]);
nPeakMap = sum(t1MapTresholded > 0, 4);

% Normalize normM0Map again because the number of components could have
% changed
normM0MapTresholded = normalizematrix(normM0MapTresholded, 4);

end

function [t1MapReduced] = dropzerocomponents(t1Map, nPeakMap)
% DROPZEROCOMPONENTS Function removes (drops) components
% (from the 4th dimension) which have 0 values for every voxel.
%
%    [t1MapReduced] = dropzerocomponents(t1Map, nPeakMap)
%
%    INPUT:
%
%	 T1map        - 4D (x, y, z, t1t) parametric map of estimated T1
%                   values.
%    nPeakMap     - 3D (x, y, z) parametric map of the number of
%                   components.
%
%    OUTPUT:
%
%    t1MapReduced - 4D (x, y, z, m) parametric map of estimated T1 values,
%                   with all-zero components removed. Parameter m is equal
%                   to themaximum number of peaks.

% Get t1Map size.
[xSize, ySize, zSize, tSize] = size(t1Map);

% Mae output array.
maxPeaks = max(nPeakMap, [], 1:3);
t1MapReduced = zeros([size(nPeakMap), maxPeaks]);

% Rehsape to long.
t1MapLong = reshape(t1Map, [xSize * ySize * zSize, tSize]);
t1MapAggLong = reshape(t1MapReduced, [xSize * ySize * zSize, maxPeaks]);

nVoxel = size(t1MapLong, 1);
for iVoxel = 1:nVoxel
    
    t1Vector = t1MapLong(iVoxel, t1MapLong(iVoxel, :) > 0);
    padSize = maxPeaks - length(t1Vector);
    t1Vector = padarray(t1Vector, [0 padSize], 0, 'post');
    t1MapAggLong(iVoxel, :) = t1Vector;
    
end

% Reshape back
t1MapReduced = reshape(t1MapAggLong, [xSize ySize zSize maxPeaks]);

end
