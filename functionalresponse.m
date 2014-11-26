function I = functionalresponse(funtype, outtype, varargin)
%FUNCTIONALRESPONSE Returns common predator-prey functional responses
%
% I = functionalresponse(funtype)
% I = functionalresponse(funtype, outtype)
% I = functionalresponse(funtype, 'inline', 'v1', 'v2', ...)
%
% Returns a function for the total ingestion of prey B by predator P
% (concentration/time).  The return value can be either a symbolic
% expression or inline function.
%
% Input variables:
%
%   funtype:    name of functional reponse.
%               'lv':           Lotka-Volterra, I = aBP
%               'lvforage':     Lotka-Volterra with foraging arena, I = aVP
%               'type2':        Holling Disk Type 2, I = aBP/(1 + haB)
%               'type2forage':  Holling Disk Type 2 with foraging arena, 
%                               I = aVP/(1 + haV)
%
%   outtype:    'sym':      returns a symbolic expression (default)
%               'inline':   returns a vectorized inline function version of
%                           the symbolic expression
%
%   v#:         Specific input arguments to be used in the inline
%               expression (only if outtype is 'inline').  These can be
%               included in order to force the input arguments to be in a
%               specific order.
%
% Output variables:
%
%   I:          1 x 1 inline function, n x 1 cell array of inline
%               functions, or n x 1 symbolic array, in terms of one or more
%               of the following variables.  Units are in terms of biomass
%               M, area or volume A, and time T.
%               B:  biomass of prey group (M/A)
%               P:  biomass of predator group (M/A)
%               a:  search parameter (A M^-1 T^-1))
%               h:  handling time (T)
%               v:  vulnerability exchange rate (T^-1)

% Copyright 2009 Kelly Kearney

%--------------------------
% Parse input
%--------------------------

if nargin < 2
    converttoinline = false;
else
    if strcmp(outtype, 'sym')
        converttoinline = false;
    elseif strcmp(outtype, 'inline')
        converttoinline = true;
    else
        error('Unknown output type');
    end
end

if ~converttoinline & nargin > 2
    warning('Named input arguments only applicable with inline function output');
end
  
%--------------------------
% Calculate symbolic 
% expression for functional
% reponse
%--------------------------

% Define the symbolic variables used in the functional responses

v = sym('v');
B = sym('B');
V = sym('V');
a = sym('a');
P = sym('P');
h = sym('h');

% Qfun is the functional response in simplest terms

switch funtype
    case 'lv'
        Qfun = a*B*P;
    case 'lvforage'
        Qfun = a*V*P;
    case 'type2'
        Qfun = (a*B*P)/(1 + h*a*B);
    case 'type2forage'
        Qfun = (a*V*P)/(1 + h*a*V);
    otherwise
        error('Unknown functional response type');
end

% Find variables used in functional response equation 

temp = textscan(findsym(Qfun), '%s', 'delimiter', ',');
qvars = temp{1};

% For foraging arena functional repsonses, solve for V (vulnerable prey) by
% assuming dV/dt = 0, and substitute back to get new expression for
% ingestion of prey

if ismember('V', qvars)
    dVdt = v*(B - 2*V) - Qfun;
    Vfun = solve(dVdt, V);
    I = subs(Qfun, V, Vfun);
else
    I = Qfun;
end

%--------------------------
% Convert to inline 
% function if necessary
%--------------------------

if converttoinline
    if length(I) == 1
        I = inline(vectorize(I), varargin{:});
    else
        Inew = cell(size(I));
        for ii = 1:length(I)
            Inew{ii} = inline(vectorize(I(ii)));
        end
        I = Inew;
    end
end

