function Param = ecosystem(Ewein, filename, tlim, dt, varargin)
%ECOSYSTEM Integrates Ecosim-like system of ODEs
%
% Param = ecosystem(Ewein, filename, tlim, dt, varargin)
%
% Input variables:
%
%   Ewein:      1 x 1 structure, Ewe input structure
%
%   filename:   string, name for output file
%
%   tlim:       1 x 2 array, simulation time limits.  Units correspond to
%               the Ewein units. 
%
%   dt:         1 x 1 array, integration time step
%
%   binit:      ngroup x 1 array, initial biomass values. If not included,
%               runs will start from the Ecopath-balanced biomass values.
%
% Optional input variables (passed as parameter/value pairs):
%
%   switching:  ngroup x 1 array, switching parameter for each predator
%               (0-2). 0 = no switching, < 1 = switches only when prey is
%               very rare, > 1 switches quickly (no units).  You can also
%               provide a scalar value to assign the same switching
%               paramter to all predators. [1]
%
%   funcres:    string indication type of functional response used for
%               predator-prey interactions.  See ecosimfeed.m for options.
%               [type2]
%
%   ppforce:    primary production forcing.  First column is time, in the
%               same units as tlim.  The remaining column lists production
%               values for all primary producers in the model, in units
%               corresponding to the Ewe input (mass area-or-volume^-1
%               time^-1). If empty, no forcing will be provided; primary
%               production will be set to balance any loss to prducer
%               groups (i.e. steady state). []
%
%   ingforce:   ingestion flux forcing.  First row hold prey index, second
%               row the predator index.  First column is time.  The
%               remaining columns list the ingestion flux from prey to
%               predator (mass area-or-volume^-1 time^-1). []


% Parameter setup 

Param = ecosystemodesetup(Ewein, varargin{:});

% Time values

tsim = tlim(1):dt:tlim(end);

nt = length(tsim);

%----------------------------
% Create output file
%----------------------------

% Start file

[blah, blah, ext] = fileparts(filename);
if isempty(ext)
    filename = [filename '.nc'];
end

mode = bitor(nc_clobber_mode, nc_64bit_offset_mode);
nc_create_empty(filename, mode);

% File description

nc_attput(filename, nc_global, 'title', 'ecosystem simulation output' );
nc_addhist(filename, 'Ecosystem simulation started');

% Add dimensions

nc_add_dimension(filename, 'time', nt);
nc_add_dimension(filename, 'group', Param.ngroup);
nss = Param.ngroup + Param.ngear + 1;
nc_add_dimension(filename, 'sourcesink', nss);

% Coordinate variables

Var.Name = 'time';
Var.Dimension = {'time'};
nc_addvar(filename, Var);
nc_varput(filename, Var.Name, tsim);

% Data variables

Var.Name = 'biomass';
Var.Dimension = {'time', 'group'};
nc_addvar(filename, Var);
nc_varput(filename, Var.Name, Param.b, [0 0], [1 Param.ngroup]); % Add starting value

if Param.diagflag
    flux = {'ingest', 'egest', 'excrete', 'fish', 'nonpred','pp'};
    nfx = length(flux);
    for ifx = 1:nfx
        Var.Name = flux{ifx};
        Var.Dimension = {'time', 'sourcesink', 'sourcesink'};
        nc_addvar(filename, Var);
    end
end
    
%----------------------------
% Simulation
%----------------------------

if isempty(Param.binit)
    b = Param.b;
else
    b = Param.binit;
end

% Loop over time

fprintf('Running ecosystem:           \n');
erasestr = repmat('\b', 1, 11);

for it = 1:(nt-1)
    
    if ~mod(it,10) 
        fprintf([erasestr sprintf('%10.2f\n', tsim(it))]);
    end
    
    % Try fixed step solver
    
    [tout, newb, db, Flux] = odewrap(@ode4,  @ecosystemode, tsim(it)+[0 dt], b, [], Param);           
    newb = endonly(newb);

    % May add in variable step on later
    
    isbad = isnan(newb) | isinf(newb) | newb < 0;
    
    if any(isbad(:))
        [tout, newb, db, Flux] = odewrap(@ode45,  @ecosystemode, tsim(it)+[0 dt], b, [], Param);           
        newb = endonly(newb);
        
        isbad = isnan(newb) | isinf(newb) | newb < 0;
        
        if any(isbad(:))
            error('NaN, Inf, or negative');
        end
    end
    
    % Write to file (fluxes correspond to beginning of time step, biomass
    % to end of time step) 

    nc_varput(filename, 'biomass', b, [it 0], [1 Param.ngroup]);

    if Param.diagflag
        for ifx = 1:nfx
            data = permute(Flux(2).(flux{ifx}), [3 1 2]);
            nc_varput(filename, flux{ifx}, data, [it-1 0 0], [1 nss nss]);
        end
    end
    
    b = newb;
    
end
fprintf('\n');

