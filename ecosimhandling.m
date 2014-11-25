function h = ecosimhandling(qb, qbmaxqb0);
%ECOSIMHANDLING Calculates handling time 
%
% h = ecosimhandling(qb, qbmaxqb0);
%
% Calculates the handling time for each predator based on Q/B and
% (Q/Bmax)/(Q/B0)
%
% Input variables:
%
%   qb:         ngroup x 1 array, consumption per unit biomass each
%               predator (time^-1) 
%
%   qbmaxqb0:   ngroup x 1 array, ratio of maximum possible consumption per
%               unit biomass compared to the mass-balance amount (no unit) 
%
% Output variables:
%
%   h:          ngroup x 1 array, handling time (time)

% Copyright 2008 Kelly Kearney

if ~isvector(qb) | ~isvector(qbmaxqb0) | length(qb)~=length(qbmaxqb0)
    error('qb and qbmaxqb0 must be vectors of the same length');
end

qb = qb(:);
qbmaxqb0 = qbmaxqb0(:);

h = (1./(qb .* qbmaxqb0));