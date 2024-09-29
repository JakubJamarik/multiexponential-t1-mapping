function [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1map(v, ti, method, varargin)
% CALCULATET1MAP Retunrs the T1 map, M0 map, normalized M0 map, and a map
% of number of components of the input 4D IR MRI image, based on the method
% selected.
%
% [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1map(v, ti, method, varargin)
%
%    INPUT:
%    v             - 4D (x, y, z, TI) matrix of IR MRI image.
% 	 ti            - Vector of TI times in ms.
%    method        - Char indicating which method should be used for T1
%                    estiamtion.
%
%                    'ilt' - The Inverse Laplace Transform method. See
%                             calculatet1mapilt for details and inputs.
%
%                    'mul' - Bounded optimization with MULtiple starting
%                             points method. See calculatet1mapmul for
%                             details and inputs.
%                    'tom' - Bounded optimization by Tomer et al. 2022.
%
%    OUTPUT:
%
%	 T1map         - 4D (x, y, z, t1t) parametric map of estimated T1
%                      values.
%    nComponentMap - 4D (x, y, z, n) parametric map indicating number
%                      of components per voxel.
%    m0Map         - 4D (x, y, z, t1t) parametric map of estimated M0
%                      values.
%    normM0Map     - 4D (x, y, z, t1t) parametric map of estimated M0
%                      values, normalized across the 4th dim.

switch method
    case 'ilt'
        %
        % Input validation
        %
        
        if length(varargin) < 1
            error('ILT requires at least 1 parameter: MIN_PEAK_HEIGHT.');
        end
        
        % Set mandatory variables
        MIN_PEAK_HEIGHT = varargin{1};
        
        %
        % Defaults
        %
        
        t1t = linspace(50, 5000, 100);
        model = 'absoluteValueOfSum';
        Defaults = {t1t, model};
        Defaults(1:(length(varargin) - 1)) = varargin(2:end);
        [t1t, model] = Defaults{:};
        m0Map = [];
        normM0Map = [];
        
        %
        % Main
        %
        
        [t1Map, nComponentMap] = calculatet1mapilt(v, ti, MIN_PEAK_HEIGHT, t1t, model);
        
    case 'mul'
        %
        % Input validation
        %
        
        if length(varargin) < 2
            error('MUL requires at least 2 parameters: model and nComponent.');
        end
        
        % Set mandatory variables
        model = varargin{1};
        nComponent = varargin{2};
        
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
        Defaults(1:(length(varargin) - 2)) = varargin(3:end);
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
        % Main
        %
        
        [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1mapmul(v, ti, model, nComponent, x0, lb, ub, options, numStartingPoints, doMultiStartInParallel, startingPointMethod);
        
        
    case 'tom'
        %
        % Input validation
        %
        
        % No mandatory parameters
        
        %
        % Defaults
        %
        
        nComponentMax = floor(size(v,4)/4); % default maximum number of components per voxel
        tolFun = 0.1;
        tolX = 0.1;
        diffMaxChange = 0.1;
        diffMinChange = 0.0001;
        maskingThreshold = 100;
        
        Defaults = {...
            nComponentMax, ...
            tolFun, ...
            tolX, ...
            diffMaxChange, ...
            diffMinChange, ...
            maskingThreshold};
        Defaults(1:length(varargin)) = varargin;
        [nComponentMax, tolFun, tolX, diffMaxChange, diffMinChange, maskingThreshold] = Defaults{:};
        
        %
        % Main
        %
        
        [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1maptom(v, ti, nComponentMax, tolFun, tolX, diffMaxChange, diffMinChange, maskingThreshold);

        
    otherwise
        error('Unknown method!');
end

end
