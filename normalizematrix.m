function [Mnormalized] = normalizematrix(M, dim)
% NORMALIZEMATRIX Normalizes all values in matrix M across the dimension
% dim to 0-1 range.
%
%    [Mnormalized] = normalizematrix(M, dim)
%
%    INPUT:
%
%    M           - Matrix to be normalized.
%    dim         - Scalar indicating whihc dimension to normalize across.
%
%    OUTPUT:
%
%    Mnormalized - Normalized matrix M across the dimension dim.

Mnormalized = M ./ sum(M, dim);
Mnormalized(isnan(Mnormalized)) = 0;
end