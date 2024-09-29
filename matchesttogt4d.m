function [estMapMatched] = matchesttogt4d(estMap, gtMap)
% MATCHESTTOGT4D Function rearanges the map of parameter estimates to best
% match the ground-truth map, across the fourth dimension. This is a
% wrapper function for matchesttogt.
%
%    [estMapMatched] = matchesttogt4d(estMap, gtMap)
%
%    INPUT:
%
%    estMap        - 4D (x, y, z, t) parametric map of estimated values.
%    gtMap         - 4D (x, y, z, t) parametric map of ground-truth values.
%
%    OUTPUT:
%
%    estMapMatched - 4D (x, y, z, t) parametric map of estimated values,
%                    matched with ground truth. 

% Reshape from 4D to 2D
[xSize, ySize, zSize, tSize] = size(estMap);
nVoxel = xSize * ySize * zSize;
estMapLong = reshape(estMap, [nVoxel tSize]);
gtMapLong = reshape(gtMap, [nVoxel tSize]);
estMapMatchedLong = zeros([nVoxel tSize]);

% Get matched estimate vector for each voxel
for iVoxel = 1:nVoxel
    
    iEstMatched = matchesttogt(estMapLong(iVoxel, :), gtMapLong(iVoxel, :));
    iEstMatched = iEstMatched';
    estMapMatchedLong(iVoxel, :) = iEstMatched;
    
end

% Reshape back to original form
estMapMatched = reshape(estMapMatchedLong, [xSize ySize zSize tSize]);

end

function [estMatched] = matchesttogt(est, gt, varargin)
% MATCHESTTOGT Function rearanges the input parameter estimate vector to best match the ground-truth vector.
% Vectors are matched based on the desired method.
%
%    [estMatched] = matchesttogt(est, gt, method)
%
%    INPUT:
%
%    est        - Vector of estimated aprameters.
%    gt         - Ground-truth vector.
%    method     - Method for est and gt Matching. The supported methods
%                 are: 'sumOfRelativeErrors', 'corelation', and 'euclidean'
%                 . Each method represents a distance metric which is used 
%                 to find the closest permutation of the est and gt vector.
%  
%    OUTPUT:
%
%    estMatched - Vector of estimated parameters matched to ground-truth
%                 vector.


%
% Defaults
%

method = "sumOfRelativeErrors";

% Overwrite defautls if supplied in varargin
Defaults = {method};
Defaults(1:length(varargin)) = varargin;
[method] = Defaults{:};

% 
% Validate input
%

if isrow(gt)
    gt = gt';
end
if isrow(est)
    est = est';
end

%
% Main
%

% Compute all permutations
estPerm = perms(est);

% Compute distance between all permutations of est and gt
switch method
    case 'sumOfRelativeErrors'
        % Compute the distance using sum of relative errors
        gtDist = sum(abs(estPerm - gt')./max(realmin, gt'), 2);
        [~, minInd] = min(gtDist);
    case 'corelation'
        % Compute distance using correlation
        nPerm = size(estPerm, 1);
        gtCorr = zeros(nPerm, 1);
        for iPerm = 1:nPerm
            gtCorr(iPerm, 1) = sum(abs(estPerm(iPerm, :)' - gt) ./ gt);
        end
        [~, minInd] = min(gtCorr);
    case 'euclidean'
        % Compute distance usngi the euclidean distance
        gtDist = pdist2(estPerm, gt', 'euclidean');
        [~, minInd] = min(gtDist);
    otherwise
        error('Method unknown!!!')
end

% Return the permutation with the lowest distance to gt
estMatched = estPerm(minInd, :);
estMatched = estMatched';

end
