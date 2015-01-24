function ing = ecosimfeed(fr, b, a, k, v, h, sw);
%ECOSIMFEED Calculate ingestion using Ecosim-based functional responses
%
% ing = ecosimfeed(fr, b, a, k, v, h, sw);
%
% Input variables:
%
%   fr:     functional response to use:
%           'lv':           Lotka-Volterra, I = aBP
%           'lvforage':     Lotka-Volterra with foraging arena, I = aVP
%           'type2':        Holling Disk Type 2, I = aBP/(1 + haB)
%           'type2forage':  Holling Disk Type 2 with foraging arena, 
%                           I = aVP/(1 + haV)
%           
%           fr can also be any inline function of the form
%           ifun(B,P,a,h,v) that returns a ngroup x ngroup array.
%
%   b:      ngroup x 1, biomass of each group (concentration)
%
%   a:      ngroup x ngroup, non-switching search parameter for each
%           predator-prey relationship (time^-1 concentration^-1).  Rows
%           correspond to prey groups and columns to predators.
%
%   k:      ngroup x ngroup, scaling factor used to convert from
%           non-switching search parameter to switching one (no units) 
%
%   v:      ngroup x ngroup, vulnerability exchange rate, (time^-1)
%
%   h:      ngroup x 1, handling time for each predator (time)
%
%   sw:     ngroup x 1, switching parameter for each predator (no unit,
%           0-2). 0 = no switching, < 1 = switches only when prey is very
%           rare, > 1 switches quickly (no units)
%
% Output variables:
%
%   ing:    ngroup x ngroup, ingestion of each prey by each predator
%           (concentration/time)

% Copyright 2008 Kelly Kearney

%-----------------------------
% Setup
%-----------------------------

% Check input

if ~isvector(b) ||  ~isvector(sw) || ~isequal(length(b),length(sw))
    error('b, h, and sw  must be vectors of the same length');
end
b = b(:);
sw = reshape(sw,1,[]);
nb = length(b);

if ~isequal([nb nb], size(v), size(a), size(k))
    error('v, a, and k must be square matrices of the same length as b');
end

if ~(isvector(h) || isequal(size(h), [nb nb]))
    error('h must be vector or square matrix of the same length as b');
end

% Reshape and resize

[bj, bi] = meshgrid(b);
if isvector(h)
    h = reshape(h,1,[]);
    h = repmat(h, nb, 1);
end
sw = repmat(sw, nb, 1);

%-----------------------------
% Ingestion calculation
%-----------------------------

% Switching

asw = bsxfun(@rdivide, k.*a.*bi.^sw, nansum(a.*bi.^sw, 1));

% Ingestion

if ischar(fr)
    switch fr
        case 'lv'
            ifun = functionalresponse('lv',       'inline', 'B', 'P', 'a', 'h', 'v');
        case 'lvforage'
            ifun = functionalresponse('lvforage', 'inline', 'B', 'P', 'a', 'h', 'v');
        case 'type2'
            ifun = functionalresponse('type2',    'inline', 'B', 'P', 'a', 'h', 'v');
        case 'type2forage'
            ifun = functionalresponse('type2',    'inline', 'B', 'P', 'a', 'h', 'v');
            ifun = ifun{1};
        otherwise
            error('Unknown functional response type');
    end
elseif strcmp(class(fr), 'inline')
    ifun = fr;
end
ing = ifun(bi, bj, asw, h, v);

% Two things can result in a NaN result: having an infinite handling time
% (from Q/B = 0 for nonconsumers) or when all prey groups of a given
% predator have 0 biomass.  In both cases, ingestion should actually be 0.

isnotconsumer = all(isinf(h) & a == 0);

bipreyonly = zeros(size(bi));
bipreyonly(a > 0) = bi(a > 0);
allpreygone = sum(bipreyonly) == 0 & ~isnotconsumer;

ing(:, isnotconsumer) = 0;
ing(:, allpreygone) = 0;    

% NOTE: The above statement only holds when using a type 2 functional
% response.  Here I set all ingestion to 0 if a = 0, since that implies no
% predation between the groups. 

nopred = a == 0;
ing(nopred) = 0;


