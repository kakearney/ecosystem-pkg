function [dbdt, Flux] = ecosystemode(t, b, A)
%ECOSYSTEMODE ODE equation for ecosystem.m
%
% This is the ODE associated with ecosystem.m, and is designed to be
% called directly by that function.

ng = A.ngroup;
ngr = A.ngear;

%--------------------------
% Calculate fluxes
%--------------------------

% All grids are square with ng + ngr + 1 elements per side, representing
% the functional groups from the Ecopath model, plus fishing gear types,
% and one extra for the mysterious source of primary production and sink
% for excretion/respiration.

[ingest, egest, excrete, fish, nonpred, pp] = deal(zeros(ng+ngr+1));
fgidx = 1:ng;

% Predation

ingest(fgidx,fgidx) = ecosimfeed(A.funcresfun, b, A.a, A.k, A.v, A.h, A.sw);
if ~isempty(A.ingforce)
    ingest = forceingest(ingest, A.ingforce, t);
end

% Egestion

detidx = find(A.pp == 2);
egest(fgidx,fgidx) = egestion(ingest(fgidx,fgidx), A.gs, A.df, detidx);

% Excretion/respiration

excrete = excretion(ingest, A.ge, A.gs);

% Fishing mortality

fish(fgidx,1:end-1) = fishing(b, A.fish);

% Other mortality, non-predatory and non-fishing related

nonpred(fgidx,fgidx) = otherloss(b, A.m0, A.df, detidx);

% Primary production

pp = primprod(A.pp, b, t, A.r, A.hpp, A.ppforce, A.pbforce, size(pp));
% 
% if ~isempty(A.ppforce)
%     pp = forcepp(pp, A.pp, A.ppforce, t);
%     if any(isnan(pp(:)))
%         pp = primprod(A.pp, ingest, egest, excrete, fish, nonpred);
%     end
% elseif ~isempty(A.pbforce)
%     pp = forcepb(A.pp, b, A.pp, A.pbforce, t)
%     if any(isnan(pp(:)))
%         pp = primprod(A.pp, ingest, egest, excrete, fish, nonpred);
%         (pptype, b, r, h, sz)
%     end
% else
% %     pp2 = primprod2(A.pp, ingest, egest, excrete, fish, nonpred);
%     pp = primprod(A.pp, b, A.r, A.hpp, size(ingest));
% end
%         

%--------------------------
% Total rate of change
%--------------------------

fluxtot = ingest + egest + excrete + fish + nonpred + pp;
fluxin = sum(fluxtot,1)';
fluxout = sum(fluxtot,2);

dbdt = fluxin(1:ng) - fluxout(1:ng);

Flux = struct(...
    'ingest', ingest, ...
    'egest', egest, ...
    'excrete', excrete, ...
    'fish', fish, ...
    'nonpred', nonpred, ...
    'pp', pp);

%**************************************************************************

%--------------------------
% Egestion
%--------------------------

function egest = egestion(ing, gs, df, detidx)

ng = size(ing,1);

suming = sum(ing, 1);

egest = zeros(ng);
for ig = 1:ng
    egest(ig, detidx) = suming(ig) * gs(ig) * df(ig,:);
end

%--------------------------
% Excretion
%--------------------------

function excrete = excretion(ing, ge, gs)
suming = sum(ing,1);
excrete = zeros(size(ing));
ng = length(ge);
excrete(1:ng,end) = suming(1:ng)' .* (1 - ge - gs);

%--------------------------
% Fishing
%--------------------------

function fish = fishing(b, frate)

ng = length(b);
ngr = size(frate,2);

fish = zeros(ng, ng + ngr);
fish(:,ng+1:end) = bsxfun(@times, b, frate);

%--------------------------
% Other mortality
%--------------------------

function nonpred = otherloss(b, m0, df, detidx)

for ig = 1:length(b)
    nonpred(ig,detidx) = b(ig) * m0(ig) * df(ig,:);
end

%--------------------------
% Primary production
%--------------------------

function pp = primprod(pptype, b, t, r, h, ppforce, pbforce, sz)

warnstat = warning('off', 'MATLAB:interp1:NaNinY');

[pp1, pp2, pp3] = deal(zeros(sz));

% Forced primary production

if isempty(ppforce)
    pp1(end, pptype == 1) = NaN;
else
    ppval = interp1(ppforce(:,1), ppforce(:,2:end), t);
    pp1(end, pptype == 1) = ppval;
end

% Forced P/B

if isempty(pbforce)
    pp2(end, pptype == 1) = NaN;
else
    pbval = interp1(pbforce(:,1), pbforce(:,2:end), t);
    ppval = pbval' .* b(pptype == 1);
    pp2(end, pptype == 1) = ppval;
end

% Unforced

pptemp = r .* b ./ (1 + h .* b);
pp3(end,pptype == 1) = pptemp(pptype == 1);

% Combine

usepp = ~isnan(pp1);
usepb = ~isnan(pp2) & isnan(pp1); % pp force takes priority
pp = pp3;
pp(usepp) = pp1(usepp);
pp(usepb) = pp2(usepb);

warning(warnstat);

% % Primary production (hack to force steady state)
% 
% function pp = primprod2(pptype, ingest, egest, excrete, fish, nonpred)
% 
% allout = sum(ingest + egest + fish + nonpred + excrete, 2);
% pp = zeros(size(ingest));
% pp(end,pptype == 1) = allout(pptype==1);
% 
% % Primary production (Ecosim formula)
% 
% function pp = primprod(pptype, b, r, h, sz)
% 
% pptemp = r .* b ./ (1 + h .* b);
% pp = zeros(sz);
% pp(end,pptype == 1) = pptemp(pptype == 1);
% 
% 
% % Forced primary production
% 
% function pp = forcepp(pp, pptype, ppforce, t)
% ppval = interp1(ppforce(:,1), ppforce(:,2:end), t);
% pp(end, pptype == 1) = ppval;
% 
% % Force primary production per unit biomass
% 
% function pp = forcepb(pp, b, pptype, pbforce, t)
% pbval = interp1(pbforce(:,1), pbforce(:,2:end), t);
% ppval = pbval' .* b(pptype == 1);
% pp(end, pptype == 1) = ppval;

%--------------------------
% Forced ingestion
%--------------------------

function ing = forceingest(ing, ingforce, t)

t2 = ingforce(3:end,1);
val = ingforce(3:end,2:end);
preyidx = ingforce(1,2:end);
predidx = ingforce(2,2:end);

val = interp1(t2, val, t);
for iv = 1:length(val)
    if ~isnan(val(iv))
        ing(preyidx(iv), predidx(iv)) = val(iv);
    end
end




