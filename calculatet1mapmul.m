function [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1mapmul(v, ti, model, nComponent, varargin)
% CALCULATET1MAPMUL Retunrs the T1 map, M0 map, normalized M0 map, and a map
% of number of components of the input 4D IR MRI image, estiamted by the
% Bounded optimization with MULtiple starting points (MUL).
%
% [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1mapmul(v, ti, model, nComponents, x0, lb, ub, options, numStartingPoints, doMultiStartInParallel, startingPointMethod)
%
%    INPUT:
%
%    v                      - 4D (x, y, z, TI) matrix of IR MRI image.
% 	 ti                     - Vector of TI times in ms.
%    model                  - Char indicating which model to use. Default
%                             is 'absoluteValueOfSum'.
%
%                             'sumOfAbsoluteValues' - The signal is modeled
%                                                     as a sum of absolute
%                                                     values of individual
%                                                     exponential
%                                                     functions.
%                              'absoluteValueOfSum' - The signal is modeled
%                                                     as an absolute value
%                                                     of a sum of
%                                                     individual
%                                                     exponential
%                                                     functions.
%    nComponents            - Scalar indicating number of components per
%                             voxel.
%    x0                     - Vector with initial values of M0 and T1, for
%                             optimization.
%    lb                     - Vector with lower bounds of M0 and T1, for
%                             optimization.
%    ub                     - Vector with upper bounds of M0 and T1, for
%                             optimization.
%    options                - A struct with optimization options as given
%                             by
%                             createOptimProblem. (see lsqcurvefit options)
%    numStartingPoints      - Scalar with number of starting points for
%                             repeated
%                             initialization of optimization.
%    doMultiStartInParallel - Logical indicating wether to run
%                             initialization of optimization in parallel.
%    startingPointMethod    - String indicating the method to for x0
%                             computation. Is computed only if x0 = [].
%                             Default is 'random'.
%
%                             'random'  - Starting point chosen at random.
%                             'tomlike' - Starting point chosen as in TOM.
%
%    OUTPUT:
%
%	 T1map           - 4D (x, y, z, t1t) parametric map of estimated T1
%                      values.
%    nComponentMap   - 4D (x, y, z, n) parametric map indicating number
%                      of components per voxel.
%    m0Map           - 4D (x, y, z, t1t) parametric map of estimated M0
%                      values.
%    normM0Map       - 4D (x, y, z, t1t) parametric map of estimated M0
%                      values, normalized across the 4th dim.

%
% Defaults
%

% Set values to [] if they will be estiamted later per voxel
x0 = [];
lb = [];
ub = [];
options = optimoptions(...
    'lsqcurvefit', ...
    'Diagnostics', 'off', ...
    'Display', 'off',...
    'FunctionTolerance', 1e-12,...
    'MaxFunctionEvaluations', 8*1000, ...
    'MaxIterations', 40000, ...
    'StepTolerance', 1e-12);

% Set defaults
Defaults = {...
    x0, ...
    lb, ...
    ub, ...
    options, ...
    10000, ...  % numStartingPoints
    false, ...  % doMultiStartInParallel
    'random' ...% startingPointMethod
    };
Defaults(1:length(varargin)) = varargin;
[...
    x0, ...
    lb, ...
    ub, ...
    options, ...
    numStartingPoints, ...
    doMultiStartInParallel, ...
    startingPointMethod, ...
    ] = Defaults{:};

%
% Input validation
%

if isrow(ti)
    ti = ti';
end

%
% Main
%

% Find size of input volume.
[xSize, ySize, zSize, tSize] = size(v);

% Reshape input volume into 2D table.
vLong = reshape(v, [xSize * ySize * zSize, tSize]);
t1MapLong = zeros(size(vLong, 1), nComponent);
m0MapLong = zeros(size(vLong, 1), nComponent);

% Compute T1 and M0 maps
nVoxel = size(vLong, 1);
parfor iVoxel = 1:nVoxel
    
    % Progress counter
    if iVoxel == 1
        disp('Fit started ')
    elseif mod(iVoxel, round(nVoxel * 0.1)) == 0
        disp('. ')
    elseif iVoxel == nVoxel
        disp('100% done.')
    end
    
    % Estimate T1 and M0 for i-th voxel
    signal = vLong(iVoxel, :)';
    if sum(signal) == 0
        disp('Skipping masked voxel.');
        m0MapLong(iVoxel, :) = 0;
        t1MapLong(iVoxel, :) = 0;
    else
        % Estimate T1 and M0 for i-th voxel
        fitTable = calculatet1mul(signal, ti, model, nComponent, x0, lb, ub, options, numStartingPoints, doMultiStartInParallel, startingPointMethod);
        
        % Reshape output
        fitTable = sortrows(fitTable, 2);  % sort based on 2. column - T1
        m0MapLong(iVoxel, :) = fitTable(:, 1)';  % first column is M0
        t1MapLong(iVoxel, :) = fitTable(:, 2)';  % second column is T1
    end
    
end

% Reshape back to 4D
t1Map = reshape(t1MapLong, [xSize ySize zSize nComponent]);
m0Map = reshape(m0MapLong, [xSize ySize zSize nComponent]);

% Compute normalized M0 map
normM0Map = normalizematrix(m0Map, 4);

% Create nComponentMap and set it to nComponent
nComponentMap = zeros([xSize, ySize, zSize]);
nComponentMap(:, :, :) = nComponent;

end

function [coefs] = calculatet1mul(signal, ti, model, nComponents, varargin)
% CALCULATET1MUL Estimates values o T1 and M0 for a multiexponential T1
% deacy signal, depending on TI.
%
% [coefs] = calculatet1mul(signal, ti, model, nComponents, x0, lb, ub, options, numStartingPoints, doMultiStartInParallel, startingPointMethod)
%
%    INPUT:
%
%    signal                 - Vector of IR signal dependent on TI.
%    ti                     - Vector of TI times in ms.
%    model                  - Char indicating which model to use. Default
%                             is 'absoluteValueOfSum'.
%
%                             'sumOfAbsoluteValues' - The signal is modeled
%                                                     as a sum of absolute
%                                                     values of individual
%                                                     exponential
%                                                     functions.
%                              'absoluteValueOfSum' - The signal is modeled
%                                                     as an absolute value
%                                                     of a sum of
%                                                     individual
%                                                     exponential
%                                                     functions.
%    nComponents            - Scalar indicating number of components per
%                             voxel.
%    x0                     - Vector with initial values of M0 and T1, for
%                             optimization.
%    lb                     - Vector with lower bounds of M0 and T1, for
%                             optimization.
%    ub                     - Vector with upper bounds of M0 and T1, for
%                             optimization.
%    options                - A struct with optimization options as given
%                             by
%                             createOptimProblem. (see lsqcurvefit options)
%    numStartingPoints      - Scalar with number of starting points for
%                             repeated
%                             initialization of optimization.
%    doMultiStartInParallel - Logical indicating wether to run
%                             initialization of optimization in parallel.
%    startingPointMethod    - String indicating the method to for x0
%                             computation. Is computed only if x0 = [].
%                             Default is 'random'.
%
%                             'random'  - Starting point chosen at random.
%                             'tomlike' - Starting point chosen as in TOM.
%
%    OUTPUT:
%
%    coefs                  - 2D matrix with nComponent rows and 2 columns,
%                             with M0 and T1 estimates in the first and the
%                             second column, respectively.

%
% Defaults
%

options = optimoptions(...
    'lsqcurvefit', ...
    'Diagnostics', 'off', ...
    'Display', 'off',...
    'FunctionTolerance', 1e-12,...
    'MaxFunctionEvaluations', 8*1000, ...
    'MaxIterations', 40000, ...
    'StepTolerance', 1e-12);
Defaults = {...
    [], ... % x0
    [], ... % lb
    [], ... % ub
    options, ...
    10000, ...  % numStartingPoints
    false, ...  % doMultiStartInParallel
    'random' ...% startingPointMethod
    };
Defaults(1:length(varargin)) = varargin;
[...
    x0, ...
    lb, ...
    ub, ...
    options, ...
    numStartingPoints, ...
    doMultiStartInParallel, ...
    startingPointMethod, ...
    ] = Defaults{:};

%
% Set function variables
%

% Set optimization function based on model
switch model
    case 'absoluteValueOfSum'
        objFun = @absolutevalueofsum;
    case 'sumOfAbsoluteValues'
        objFun = @sumofabsolutevalues;
    otherwise
        error('Unknown model selected!')
end

% Set x0 if none was provided
if isempty(x0)
    switch startingPointMethod
        case 'random'
            x0 = [rand(1, nComponents) * max(signal) rand(1, nComponents) * 4000 / 2];
        case 'tomlike'
            x0 = [rand(1, nComponents) * max(signal) 1200 linspace(600, 2500, nComponents - 1)];
        otherwise
            error('Unknown method specified!')
    end
end

% Set lb if none was provided
minM0 = max(signal) * 0.05; % Lowest relative weight of 5%
minT1 = 10;
if isempty(lb)    
    lb = [repelem(minM0, nComponents) repelem(minT1, nComponents)];    
end

% Set ub if none was provided
maxM0 = max(signal) * 0.95; % Maximum relative weight of 95%
maxT1 = 4000;
if isempty(ub)    
    ub = [repelem(maxM0, nComponents) repelem(maxT1, nComponents)];    
end

% Set optimization problem
problem = createOptimProblem('lsqcurvefit', ...
    'objective', objFun,...
    'xdata', ti, ...
    'ydata', signal,...
    'x0', x0, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);
ms = MultiStart('UseParallel', doMultiStartInParallel);

%
% Main
%

% Perform optimization
x = run(ms, problem, numStartingPoints);

% Reshape output
coefs(:, 1) = x(1:nComponents);
coefs(:, 2) = x(nComponents+1:nComponents*2);

end

function funValue = absolutevalueofsum(p, ti)
nComponents = round(length(p) / 2);
m0 = p(1:nComponents);
t1 = p(nComponents+1:end);

t1 = 1 ./ t1;
K = - ti * t1;
K = 1 - 2 * exp(K);
K = m0 .* K;

funValue = sum(K, 2);
funValue = abs(funValue);

end

function funValue = sumofabsolutevalues(p, ti)

nComponents = round(length(p) / 2);
m0 = p(1:nComponents);
t1 = p(nComponents+1:end);

t1 = 1 ./ t1;
K = - ti * t1;
K = 1 - 2 * exp(K);
K = abs(K);
K = m0 .* K;

funValue = sum(K, 2);

end