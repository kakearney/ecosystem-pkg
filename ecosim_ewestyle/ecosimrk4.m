function [t,y] = ecosimrk4(ystart, tstart, tend, step, derivfun)
%ECOSIMRK4 Implements Ecosim-modified Runge-Kutta
%
% [t,y] = ecosimrk4(ystart, tstart, tend, step, derivfun)
%
% This function implements a 4th order Runge-Kutta integration with a fixed
% step size, with modifications for groups where the given step size does
% not resolve the changes in biomass.  When inflow/biomass for time t is
% greater than 6, biomass(t+1) = (inflow/outlow) * biomass(t).
% 
% Input variables:
%
%   ystart:     ngroup x 1 array, starting biomass values
%
%   tstart:     scalar, starting time
%
%   tend:       scalar, ending time
%
%   step:       scalar, time step size
%
%   derivfun:   handle to the main ecosim function; this calculates dB/dt
%               for and given time and biomass value, and also returns many
%               of the intermediate variables used to calculate this value
%
% Output variables:
%
%   t:          times
%
%   y:          biomass values at each time step

% Copyright 2007 Kelly Kearney


%-----------------------------
% Setup
%-----------------------------

t = tstart:step:tend;
nt = length(t);
nstep = nt - 1;

y = zeros(length(ystart), nt);
y(:,1) = ystart;

%-----------------------------
% Loop over each time step
%-----------------------------

for istep = 1:nstep
    
    % Use Runge-Kutta integration for most

    t1 = t(istep);
    y1 = y(:,istep);
    ytemp = rkstep(t(istep), step, y(:,istep), derivfun);
    
    % Use inflow/outlow approximation for rapidly-changing groups
    
    [dbdt, fishloss, otherloss, q, qin, qout, production, detex, ...
        det1a, det1b, det2] = feval(derivfun, t1, y1, true);
    
    isdet = false(size(dbdt));
    isdet(end-length(detex):end) = true; % Detritus groups always last
                                       
    inflow = production;
    outflow = qout + (fishloss + otherloss) .* y1;
    outflow(isdet) = outflow(isdet) + detex;
    
    useapprox = production./y1 > 6;
    ytemp(useapprox) = inflow(useapprox) ./ outflow(useapprox) .* y1(useapprox);
    
    y(:,istep+1) = ytemp;
    
end

%-----------------------------
% Subfunction: Compute one
% time step using Runge-Kutta
% algorithm
%-----------------------------

function y2 = rkstep(t1, h, y1, derivfun)

k1 = h * feval(derivfun, t1,       y1);
k2 = h * feval(derivfun, t1 + h/2, y1 + k1/2);
k3 = h * feval(derivfun, t1 + h/2, y1 + k2/2);
k4 = h * feval(derivfun, t1 + h,   y1 + k3);

dy = k1/6 + k2/3 + k3/3 + k4/6;
y2 = y1 + dy;


