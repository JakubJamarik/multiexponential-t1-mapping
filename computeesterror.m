function [errorMap] = computeesterror(estMap, gtMap, varargin)
% COMPUTEESTRERROR Computes error between a parameter estiamte map
% and a ground-truth map, using one of the error metrics provided.
%
%    [errorMap] = computesterror(estMap, gtMap, method, errorConstat)
%
%    INPUT:
%
%    estMap        - 4D (x, y, z, t) matrix of parameter estimate.
%    gtMap         - 4D (x, y, z, t) matrix of parameter grpund-truth.
%    method        - String idiciating which method should be used to
%                    compute the error. The default is 'relativeerror':
%
%                    'relativeerror'    - error computed as relative error
%                                         in percent.
%                    'relativeerrorwc'  - error computed as relative error
%                                         in percent with a small postive
%                                         constant added to ground-truth.
%                    'maximummeanerror' - error computed as modified mximum
%                                         mean error.   
%
%    errorConstant - Scalar inidcating value of the constat in
%                    'relativeerrorwc' method. Default is 0.001.
%
%    OUTPUT:
%
%    errorMap - 4D (x, y, z, t) matrix of maximum mean errors.

%
% Defaults
%

method = 'relativeerror';
errorConstant = 0.001;

% Overwrite defautls if supplied in varargin
Defaults = {method};
Defaults(1:length(varargin)) = varargin;
[method] = Defaults{:};

%
% Main
%

% If estimate has less values then ground-truth pad it with zeroes
if (size(gtMap, 4) > size(estMap, 4))
    padSize = size(gtMap, 4) - size(estMap, 4);
    estMap = padarray(estMap, [0 0 0 padSize], 0, 'post');
else
    % All good.
end

% Compute desired
switch method
    case 'relativeerror'
        errorMap = abs(gtMap - estMap) ./ gtMap * 100;
    case 'relativeerrorwc'
        errorMap = abs(gtMap - estMap) ./ (errorConstant + gtMap) * 100;
    case 'maximummeanerror'
        errorMap = computemmer(estMap, gtMap);
    otherwise
        error('Unknown method!!!')
end

end

function [errorMap] = computemmer(estMap, gtMap)
% COMPUTEMMER Computes Maximum Mean Error (MMER) of an estiamte map with
% respect to a ground-truth map. The MMER is given by the follwing
% equation:
%
% MMER = 0                     , x = y = 0
% MMER = abx(x - y) / max(x, y), otherwise
%
% where x is an element of the estiamte map and y is an element of the
% ground-truth map.
%
%    [errorMap] = computemmer(estMap, gtMap)
%
%    INPUT:
%
%    estMap   - 4D (x, y, z, t) matrix of parameter estimate.
%    gtMap    - 4D (x, y, z, t) matrix of parameter grpund-truth.
%
%    OUTPUT:
%
%    errorMap - 4D (x, y, z, t) matrix of maximum mean errors.

% Get the size of the last dimension (the dimension holding the values)
lastDim = ndims(gtMap);

% Initialize the error matrix
errorMap = zeros(size(gtMap));

% Loop over all elements in all dimensions except the last one
% Use indexing with the last dimension
for idx = 1:numel(gtMap)
    % Linear index converted to subscript for all dimensions
    subscripts = cell(1, lastDim);
    [subscripts{:}] = ind2sub(size(gtMap), idx);
    
    % Extract the current ground truth and estimate values along the last dimension
    x = estMap(subscripts{:});
    y = gtMap(subscripts{:});
    
    if x == 0 && y == 0
        errorMap(subscripts{:}) = 0;
    else
        errorMap(subscripts{:}) = abs(x - y) / max(x, y);
    end
end

end
