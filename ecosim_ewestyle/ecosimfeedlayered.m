function [ingin, ingout] = ecosimfeedlayered(fr, b, a, k, v, h, sw, feedtype)
%ECOSIMFEEDLAYERED Calculate ingestion using Ecosim-based functional 
%                  responses with layers
%
% [ingin, ingout] = ecosimfeedlayered(fr, b, a, k, v, h, sw, feed)
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
%   b:      ndepth x ngroup, biomass of each group (concentration)
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
%   feed:   ngroup x ngroup array of feeding types, where 1 indicates
%           nekton preying on nekton, 2 indicated nekton preying on
%           plankton, and 3 indicates plankton preying on plankton.
%
% Output variables:
%
%   ingin:  ngroup x ngroup x ndepth, ingestion of each prey by each
%           predator at each depth
%
%   ingout: ngroup x ngroup x ndepth, predation loss for each prey due to
%           each predator at each depth.  The two output array will be
%           identical for nekton-nekton and plankton-plankton fluxes, but
%           nekton-plankton feeding flows from many layers to just one.

% Copyright 2009 Kelly Kearney

%----------------------------
% Check input
%----------------------------

[ndepth, nb] = size(b);

%----------------------------
% Nekton eat nekton
%----------------------------

bnn = b(1,:)';

ingnn = ecosimfeed(fr, bnn, a, k, v, h, sw);
ingnn(feedtype ~= 1) = 0;

ingnnspread = zeros(nb, nb, ndepth);
ingnnspread(:,:,1) = ingnn;

%----------------------------
% Nekton eat plankton
%----------------------------

% Calculate using water column sum

bnp = nansum(b, 1);                 % total over water column
btemp = b;
btemp(isnan(btemp)) = 0;
bioFrac = bsxfun(@rdivide, btemp, sum(btemp,1)); % fraction in each layer
bioFrac(isnan(bioFrac)) = 0;    % occurs if b = 0 for a group

ingnp = ecosimfeed(fr, bnp, a, k, v, h, sw);
ingnp(feedtype ~= 2) = 0;

% Distribute losses proportionately over layers

ingnpspread = bsxfun(@times, ingnp, permute(bioFrac, [2 3 1]));
% 
% ingnpspread = zeros(nb, nb, ndepth); 
% for iprey = 1:nb
%     for ipred = 1:nb
%         for idepth = 1:ndepth
%             ingnpspread(iprey,ipred,idepth) = ingnp(iprey,ipred) .* bioFrac(idepth,iprey);
%         end
%     end
% end
    
%----------------------------
% Plankton eat plankton
%----------------------------

ingpp = zeros(nb, nb, ndepth);
for idepth = 1:size(b,1)
    bpp = b(idepth,:)';
    ingpptemp = ecosimfeed(fr, bpp, a, k, v, h, sw);
    ingpptemp(feedtype ~= 3) = 0;
    ingpp(:,:,idepth) = ingpptemp;
end

% 
% for idepth = 1:size(b,1)
%     bpp = b(idepth,:)';
%     ingpp{idepth} = ecosimfeed(fr, bpp, a, k, v, h, sw);
%     ingpp{idepth}(feedtype ~= 3) = 0;
% end
% ingpp = cat(3, ingpp{:});


%----------------------------
% Combine to get total 
% ingestion
%----------------------------

ingout = ingnnspread + ingnpspread + ingpp;
ingin = ingpp;
ingin(:,:,1) = ingin(:,:,1) + ingnp + ingnn;

