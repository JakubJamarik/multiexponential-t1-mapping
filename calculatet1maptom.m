function [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1maptom(v, ti, varargin)
% CALCULATET1MAPTOM Retunrs the T1 map, M0 map, normalized M0 map, and a map
% of number of components of the input 4D IR MRI image, estiamted by the
% Bounded optimization by Tomer et al. 2022 (TOM).
%
% This is a wrapper function for the original code from Tomer, O.,
% Barazany, D., Baratz, Z., Tsarfaty, G. & Assaf, Y.
% In vivo measurements of lamination patterns in the human
% cortex. Human Brain Mapping 43, 2861–2868 (2022). The original code was
% acquired from https://github.com/omritomer/t1_layers on or before
% 2024-09-16.
%
%    [t1Map, nComponentMap, m0Map, normM0Map] = calculatet1maptom(v, ti, nComponentMax, tolFun, tolX, diffMaxChange, diffMinChange, maskingThreshold)
%
%    INPUT:
%
%    v                - 4D (x, y, z, TI) matrix of IR MRI image.
% 	 ti               - Vector of TI times in ms.
%    nComponentMax    - Scalar indicating number of components per voxel.
%    tolFun           - Scalar indicating termination tolerance on the
%                       optimized function value. Default is 0.1. For more
%                       see option *FunctionTolerance* in lsqnonlin.
%    tolX             - Scalar indicating termination tolerance on
%                       optimizaton function input.
%                       Default is 0.1. For more see option *StepTolerance*
%                       in lsqnonlin.
%    diffMaxChange    - Scalar indicating the maximum change in variables
%                       for finite difference gradients. Default is 0.1.
%                       For more see lsqnonlin.
%    diffMinChange    - Scalar indicating minimum change in variables for
%                       finite difference gradients. Default is 0.0001.
%                       For more see lsqnonlin.
%    maskingThreshold - Scalar indicating threshold for masking the data.
%                       Voxel's whose mean signal is lower than threshold
%                       are treated as having no signal and therefore no
%                       T1-components are calculated for them.
%
%    OUTPUT:
%
%	 T1map            - 4D (x, y, z, t1t) parametric map of estimated T1
%                       values.
%    nComponentMap    - 4D (x, y, z, n) parametric map indicating number
%                       of components per voxel.
%    m0Map            - 4D (x, y, z, t1t) parametric map of estimated M0
%                       values.
%    normM0Map        - 4D (x, y, z, t1t) parametric map of estimated M0
%                       values, normalized across the 4th dim.

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
% Load functions
%

%
% Main
%

[t1Map, m0Map, normM0Map, nComponentMap, ~, ~, ~, ~] = calcT1map( ...
    v, ...
    ti, ...
    'Nmax', nComponentMax, ...
    'TolFun', tolFun, ...
    'TolX', tolX, ...
    'DiffMaxChange', diffMaxChange, ...
    'DiffMinChange', diffMinChange, ...
    'MaskingThreshold', maskingThreshold);

end