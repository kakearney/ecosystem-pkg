function A = ecosystemodesetup(Ewein, varargin)

%---------------------------
% Parse and check input
%---------------------------

Ewein = ecopathinputcheck(Ewein, true);
npp = sum(Ewein.pp == 1);

p = inputParser;

p.addParamValue('switching', 1, @(x) isscalar(x) || (isvector(x) && length(x) == Ewein.ngroup));
p.addParamValue('funcres', 'type2', @(x) ismember(x, {'lv', 'lvforage', 'type2', 'type2forage'}));
p.addParamValue('ppforce', [], @(x) size(x,2) == npp+1);
p.addParamValue('pbforce', [], @(x) size(x,2) == npp+1);
p.addParamValue('ingforce', [], @(x) all(reshape(ismember(x(1:2,2:end), 1:Ewein.ngroup),[],1)));
p.addParamValue('binit', [], @(x) isvector(x) && length(x) == Ewein.ngroup);
p.addParamValue('diagflag', false, @(x) islogical(x) && isscalar(x));

p.parse(varargin{:});
Opt = p.Results;

if isscalar(Opt.switching)
    Opt.switching = ones(Ewein.ngroup,1) * Opt.switching;
end


%---------------------------
% Calculate variables
% needed for ecosystemode
%---------------------------

% Run Ecopath to get initial conditions

Ep = ecopathlite(Ewein);

% Fields taken directly

A.ngroup = Ewein.ngroup;    % # of functional groups
A.ngear = Ewein.ngear;      % # of fishing gear types
A.name = Ewein.name;        % Names of functional groups
A.b = Ep.b;                 % Biomass       
A.ge = Ep.ge;               % Growth efficiency
A.gs = Ewein.gs;            % Fraction egested         
A.df = Ewein.df;            % Detritus fate
A.m0 = Ep.otherMortRate;    % Non-predatory, non-fishing mortality 
A.pp = Ewein.pp;            % Primary producer indicator
A.ppforce = Opt.ppforce;    % Primary production forcing
A.pbforce = Opt.pbforce;    % P/B ratio forcing
A.ingforce = Opt.ingforce;  % Ingestion forcing
A.binit = Opt.binit;        % Initial biomsas
A.diagflag = Opt.diagflag;  % Flag to output diagnistic variables

% Handling time (T)

A.h = ecosimhandling(Ep.qb, Ewein.qbmaxqb0);

% Vulnerability (T^-1)

A.q0 = Ep.q0;
A.q0(:, Ewein.pp==2) = 0; % Flow to detritus is NOT consumption

A.v = ecosimvulnerability(Ewein.kv, A.q0, Ep.b);

% Switching factor (no unit)

A.sw = Opt.switching;

% Search factor (A M^-1 T^-1) and scaling factor (no unit)

[A.a, A.k] = ecosimsearch(Opt.funcres, Ep.b, A.q0, A.v, A.h, A.sw);

A.a(isnan(A.a) | isinf(A.a)) = 0;
A.k(isnan(A.k) | isinf(A.k)) = 0;

% Maximum relative P/B for primary producers

ispp = Ewein.pp == 1;

A.r = zeros(Ewein.ngroup,1);
A.r(ispp) = Ewein.maxrelpb(ispp) .* Ep.pb(ispp);

% h related to maximum net primary production (Walters et.al. 1997, pg 144)
    
A.hpp = zeros(Ewein.ngroup,1);
A.hpp(ispp) = (A.r(ispp)./Ep.pb(ispp) - 1) ./ Ep.b(ispp);

% Fishing 

A.fish = Ep.fishMortRate;
A.fish(isnan(A.fish)) = 0;

if any(Ewein.discard(:) > 0)
    warning('No discarding yet; assuming all landed');
end

% Function for feeding

A.funcresfun = functionalresponse(Opt.funcres, 'inline', 'B', 'P', 'a', 'h', 'v');


