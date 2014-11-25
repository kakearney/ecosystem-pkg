function [a, k] = ecosimsearch(fr, b, q, v, h, sw)
%ECOSIMSEARCH Calculates the search parameters for ecosim feeding
%
% [a, k] = ecosimsearch(fr, b, q, v, h, sw)
%
% Input variables:
%
%   fr: functional response to use:
%       'lv':           Lotka-Volterra, I* = aB
%       'lvforage':     Lotka-Volterra with foraging arena, I* = aV
%       'type2':        type 2, I* = aB/(1 + haB)
%       'type2forage':  type 2 with foraging arena, I* = aV/(1 + haV)
%
%   b:  ngroup x 1, mass-balanced biomass (concentration)
%
%   q:  ngroup x ngroup, mass-balanced consumption of each prey (rows) by
%       each predator (columns) (concentration/time) 
%
%   v:  ngroup x ngroup, vulnerability exchange rate, (time^-1)
%
%   h:  ngroup x 1, handling time for each predator (time)
%       OR
%       ngroup x ngroup, handling time for each predator-prey relationship
%       (time)
%
%   sw: ngroup x 1, switching parameter for each predator (no unit, 0-2).
%       0 = no switching, < 1 = switches only when prey is very rare, > 1
%       switches quickly (no units)
%
% Output variables:
%
%   a:  ngroup x ngroup, non-switching search parameter for each
%       predator-prey relationship (time^-1 concentration^-1)
%
%   k:  ngroup x ngroup, scaling factor used to convert from non-switching
%       search parameter to switching one (no units)

% Copyright 2008 Kelly Kearney

%-----------------------------
% Setup
%-----------------------------

% Check input

if ~isvector(b) || ~isvector(sw) || ~isequal(length(b),length(sw))
    error('b, and sw  must be vectors of the same length');
end
b = b(:);
sw = reshape(sw,1,[]);
nb = length(b);

if ~isequal([nb nb], size(v), size(q))
    error('v and q must be square matrices of the same length as b');
end

if ~(isvector(h) || isequal(size(h), [nb nb]))
    error('h must be vector or square matrix of the same length as b');
end

% Reshape and resize 

[bj, bi] = meshgrid(b);
sw = repmat(sw, nb, 1);
if isvector(h)
    h = reshape(h,1,[]);
    h = repmat(h, nb, 1);
end

%-----------------------------
% Calculate parameters
%-----------------------------

istar = q./bj;

%  Calculate search factor for no-switching definition

switch fr
    case 'lv'
        a =  istar./bi;
    case 'lvforage'
        a = 2.*istar.*v./(-istar.*bj+v.*bi);
    case 'type2'
        a = -istar./bi./(istar.*h-1);
    case 'type2forage'
        a = -2.*istar.*v./(-istar.^2.*h.*bj-v.*bi+istar.*h.*v.*bi+istar.*bj);
    otherwise
        error('Unknown functional response type');
end

% Calculate k and a for switching definition

absw = a .* bi.^sw;
abswsum = nansum(absw, 1);
abswsum = repmat(abswsum, nb, 1);
abswratio = absw./abswsum;
k = a./abswratio;

% asum = nansum(a .* bi.^sw, 1);
% asum = repmat(asum, nb, 1);
% asw = a./asum;
% k = a./asw;
% k = repmat(asum, nb, 1)./(a .* bi.^sw);

    