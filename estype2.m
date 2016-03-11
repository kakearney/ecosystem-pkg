function ingest = estype2(a, b, v, h, k, p)
%ESTYPE2 Holling type 2 functional response as implemented in Ecosim
%
% ingest = estype2(a, b, v, h, k, p)
%
% This function calculates total ingestion of each prey by each predator,
% using the Holling type 2 response similar to that implemented within
% Ecosim, including foraging arena dynamics and switching behavior, but not
% the seasonal forcing or mediation.
%
% Input variables:
%
%   a:      search rate, n x n, prey x pred (time^-1 concentration^-1)
%
%   b:      biomass, n x 1 (concentration = mass/area or mass/volume)
%
%   v:      vulnerability exchange rate, n x n (time^-1)
%
%   h:      handling time, n x n (time)
%
%   k:      scaling constant, n x n (no units)
%
%   p:      switch factor (0-2), n x 1, 0 = no switching, < 1 = switches
%           only when prey is very rare, > 1 switches quickly (no units)
%
% Output variables:
%
%   ingest: ingestion amount (concentration/time)

% Copyright 2007 Kelly Kearney

b = reshape(b, [], 1); % corresponds to prey
p = reshape(p, 1, []); % corresponds to predator

vulprey = bsxfun(@times, v, b) ./ (2.*v + bsxfun(@times, a, b'));

bswitch = bsxfun(@power, b, p);

aswitch = bsxfun(@rdivide, k .* a .* bswitch, nansum(a .* bswitch, 1));

ingest = aswitch .* vulprey./(1 + h .* aswitch .* vulprey);

ingest = bsxfun(@times, ingest, b');