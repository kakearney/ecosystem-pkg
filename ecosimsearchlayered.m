function [afinal, k] = ecosimsearchlayered(fr, b, q, vv, hh, sw, frac, feedtype)
%ECOSIMSEARCHLAYERED Calculates the search parameters for ecosim feeding
%                    with plankton/nekton distinctions
%
% [a, k] = ecosimsearchlayered(fr, b, q, v, h, sw, frac, feed)
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
%   b:      ngroup x 1, mass-balanced biomass (concentration)
%
%   q:      ngroup x ngroup, mass-balanced consumption of each prey (rows)
%           by each predator (columns) (concentration/time)  
%
%   v:      ngroup x ngroup, vulnerability exchange rate, (time^-1)
%
%   h:      ngroup x 1, handling time for each predator (time)
%
%   sw:     ngroup x 1, switching parameter for each predator (no unit,
%           0-2). 0 = no switching, < 1 = switches only when prey is very
%           rare, > 1 switches quickly (no units) 
%
%   frac:   ndepth x 1 array, fraction of total planktonic biomass in each
%           depth layer  
%
%   feed:   ngroup x ngroup array of feeding types, where 1 indicates
%           nekton preying on nekton, 2 indicated nekton preying on
%           plankton, and 3 indicates plankton preying on plankton.
%
% Output variables:
%
%   a:      ngroup x ngroup, non-switching search parameter for each
%           predator-prey relationship (time^-1 concentration^-1)
%
%   k:      ngroup x ngroup, scaling factor used to convert from
%           non-switching search parameter to switching one (no units) 

% Copyright 2009 Kelly Kearney

%-----------------------------
% Setup
%-----------------------------

% Check input

if ~isvector(b) || ~isvector(sw) || ~isequal(length(b),length(sw))
    error('b and sw  must be vectors of the same length');
end
b = b(:);
sw = reshape(sw,1,[]);
nb = length(b);

if ~isequal([nb nb], size(vv), size(q))
    error('v and q must be square matrices of the same length as b');
end

if ~(isvector(hh) || isequal(size(hh), [nb nb]))
    error('h must be vector or square matrix of the same length as b');
end

% Reshape and resize 

[bj, bi] = meshgrid(b);
sw = repmat(sw, nb, 1);

if isvector(hh)
    hh = reshape(hh,1,[]);
    hh = repmat(hh, nb, 1);
end

% Check depth fraction

tempfrac = roundn(frac, -2);
if (sum(tempfrac) - 1) > .001
    tempfrac = round(frac/.005)*.005;
    if (sum(tempfrac) - 1) > .001
        error('Sum of frac vector must be 1');
    end
end
frac = tempfrac;

%------------------------------
% Get a symbolic expression for 
% the search rate a, assuming
% all prey is planktonic
%------------------------------

% Basic functional response

v = sym('v');
B = sym('B');
V = sym('V');
a = sym('a');
P = sym('P');
h = sym('h');

I = functionalresponse(fr, 'sym');

% Nekton-eat-nekton and nekton-eat-plankton follow the basic response of
% one big population eating one big population, but plankton-eat-plankton
% consists of lots of independant little populations eating little
% populations.

Inn = I;
for id = 1:length(frac)
    tmp(id) = subs(I, {B, P}, {B*frac(id), P*frac(id)});
end
Ipp = sum(tmp);

% Solve for a in each ingestion equation (may be multiple solutions
% depending on biomass distribution). 

Q = sym('Q');
ann = solve(Q - Inn, a);
app = solve(Q - Ipp, a);

nsol = [length(ann) length(app)]; % length(anp)
if any(nsol > 1)
    warning('KAK:layeredSolve', 'Haven''t figured out how to choose symbolic "a" solution,\n Choose one below, or try simplifying the initial biomass distribution for plankton');
    if nsol(1) > 1
        ann = ann(1);
    end
    if nsol(2) > 1
        disp(app);
        isol = input('Enter number of solution to use: app = ');
        app = app(isol);
    end
end

%-----------------------------
% Calculate parameters
%----------------------------- 

% Convert symbolic expression to function

annfun = inline(vectorize(ann), 'Q', 'B', 'P', 'h', 'v');
appfun = inline(vectorize(app), 'Q', 'B', 'P', 'h', 'v');

ann = annfun(q, bi, bj, hh, vv);
app = appfun(q, bi, bj, hh, vv);

afinal = zeros(size(q));
afinal(feedtype == 1 | feedtype == 2) = ann(feedtype == 1 | feedtype == 2);
afinal(feedtype == 3) = app(feedtype == 3);

% Calculate k and a for switching definition

absw = afinal .* bi.^sw;
abswsum = nansum(absw, 1);
abswsum = repmat(abswsum, nb, 1);
abswratio = absw./abswsum;
k = afinal./abswratio;

