function X = nf_modelfun_ASW (varargin)

% X = nf_modelfun_ASW (T,S,P,t,tracers)
% X = nf_modelfun_ASW (info)
%
% Returns gas concentrations and isotope ratios in air saturated water (ASW) or usage information about the function for NOBLEGASFIT, see OUTPUT as described below.
%
% *** OPTION-A: function X = nf_modelfun_ASW (T,S,P,t,tracers)
%
% INPUT
% T: water temperature (deg.C)
% S: salinity (g/kg)
% P: atmospheric pressure (hPa, which is the same as millibars)
% t: date of gas equilibration (decimal calendar year, i.e., t=2000.5 corresponds to the middle of year 2000). t may be left empty (t = []) for use with tracers that are not transient (the default value will be used), but specifying a (dummy) value will speed up evaluation. 
% tracers: list of tracers for which the output should be calculated (cell string)
%
% OUTPUT:
% X: output data (vector of struct with fieldnames as given in 'tracers' at input). Gas concentrations are given in the same format asin the output of nf_atmos_gas.m . Concentration ratios are given as (for example) X.RHe: 3He/4He, X.RNe: 20Ne/22Ne, X.RAr: 36Ar/40Ar.
%
% *** OPTION-B: function X = nf_modelfun_ASW (info)
%
% INPUT:
% info: string indicating the type of information needed:
%   info = 'PARAMS': returns parameter names used by this function
%   info = 'DEFAULTS': returns default values for initial values to be used for fitting
%   info = 'RANGES': returns default parameter value ranges allowed for fitting
%   info = 'TRACERS': returns the tracers included in the model
%
% OUTPUT:
% X: struct with usage information:
%   info = 'PARAMS' parameter names (cell string)
%   info = 'DEFAULTS' returns struct with default parameter values
%   info = 'RANGES' returns struct with defaults for min/max parameter values
%   info = 'TRACERS' returns cell string with names of the tracers included in the model function
%
% *******************************************************************
% This file is part of NOBLEGASFIT.
% 
% NOBLEGASFIT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% NOBLEGASFIT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with NOBLEGASFIT.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) 2014 Matthias S. Brennwald, matthias.brennwald@eawag.ch
% Additional information about NOBLEGASFIT is available at http://homepages.eawag.ch/~brennmat/
% Please contribute if you find this software useful.
% *******************************************************************

if nargin == 1
    switch upper(varargin{1})

        case 'PARAMS'
            X = { 'T' , 'S' , 'P' , 't' };

        case 'DEFAULTS'
            X.T = 10;
            X.S = 0;
            X.P = exp (-500/8300) * 1013.25;
            X.t = 2005.5; % middle of year 2005

        case 'RANGES'
            X.T.min = 0;
            X.T.max = 100;
            X.S.min = 0;
            X.S.max = 100;
            X.P.min = 0;
            X.P.max = 1013.25;
            X.t.min = -Inf;
            X.t.max = +Inf;

        case 'TRACERS'
            X = { 'He',  'Ne' , 'Ar' , 'Kr',  'Xe' , 'N2' , 'RHe' , 'RNe' , 'RAr' , 'He_3' , 'He_4' , 'Ne_20' , 'Ne_22' , 'Ar_36' , 'Ar_40' , 'Kr_84' , 'Kr_86' , 'Xe_136' , 'N2', 'SF6' };

        otherwise
            error (sprintf('nf_modelfun_ASW: unknown usage key ''%s''.',varargin{1}));

    end
else
    T = varargin{1};
    S = varargin{2};
    P = varargin{3};
    t = varargin{4};
    if isempty(t)
        t = getfield (nf_modelfun_ASW('DEFAULTS'),'t');
    end
    tracers = varargin{5};
    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end
    tracers = strrep (tracers,'-','_'); % make sure entries like He-3 work out as He_3

    if T < 0
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_ASW: T=%g is negative. This is not physically sensible in the ASW model...',T));
    end
    
    if S < 0
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_ASW: S=%g is negative. This is not physically sensible in the ASW model...',S) );
    end

    X = repmat (NA,1,length(tracers));

    for i = 1:length(tracers)
        switch tracers{i}
        
            case 'RHe'
                X(i) = nf_atmos_gas ('He-3',T,S,P,t) / nf_atmos_gas ('He-4',T,S,P,t);
                
            case 'RNe'
                X(i) = nf_atmos_gas ('Ne-20',T,S,P,t) / nf_atmos_gas ('Ne-22',T,S,P,t);
            
            case 'RAr'
                X(i) = nf_atmos_gas ('Ar-36',T,S,P,t) / nf_atmos_gas ('Ar-40',T,S,P,t);
            
            otherwise
                X(i)  = nf_atmos_gas (tracers{i},T,S,P,t);
                
        end % switch
    end % for
end
