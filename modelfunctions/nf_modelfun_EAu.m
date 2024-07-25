function X = nf_modelfun_EAu (varargin)

% X = nf_modelfun_EAu (A,t,tracers)
% X = nf_modelfun_EAu (info)
%
% Returns gas concentrations resulting from complete dissolutin of air (unfractionated excess air), see OUTPUT as described below.
%
% *** OPTION-A: function X = nf_modelfun_EAu (A,t,tracers)
%
% INPUT
% t, tracers: see nf_modelfun_ASW
% A: amount of dry air injetected per unit mass of water, "excess air concentration" (ccSTP-air / g-water)
%
% OUTPUT:
% X: output data (vector of struct with fieldnames as given in 'tracers' at input), see also nf_modelfun_ASW.m
%
% *** OPTION-B: function X = nf_modelfun_EAu (info)
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
            X = { 'A' , 't' };

        case 'DEFAULTS'
            X.A = 1E-3; % 1E-3 ccSTP-air / g-water
            X.t = getfield (nf_modelfun_ASW('DEFAULTS'),'t'); % use same value as with nf_modelfun_ASW
            
        case 'RANGES'
            X.A.min = 0;
            X.A.max = 1;
            X.t     = getfield (nf_modelfun_ASW('RANGES'),'t'); % use same values as with nf_modelfun_ASW

        case 'TRACERS'
            X = { 'He',  'Ne' , 'Ar' , 'Kr',  'Xe' , 'N2' , 'RHe' , 'RNe' , 'RAr' , 'He_3' , 'He_4' , 'Ne_20' , 'Ne_22' , 'Ar_36' , 'Ar_40' , 'Kr_84' , 'Kr_86' , 'Xe_136' , 'N2' , 'SF6' };

        otherwise
            error (sprintf('nf_modelfun_ASW_EAu: unknown usage key ''%s''.',varargin{1}));

    end
else
    A = varargin{1};
    t = varargin{2};
    if isempty(t)
        t = getfield (nf_modelfun_EAu('DEFAULTS'),'t');
    end
    tracers = varargin{3};
    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end
    tracers = strrep (tracers,'-','_'); % make sure entries like He-3 work out as He_3

    if A < 0
        warning ('noblefit:modelfun_params_out_of_bounds',sprintf('nf_modelfun_EAu: A=%g is negative. This is not physically sensible in the unfractionated-air excess air model...',A));
    end

    X = repmat (NA,1,length(tracers));

    for i = 1:length(tracers)
        switch tracers{i}
        
            case 'RHe'
                X(i) = nf_modelfun_EAu (A,t,'He-3') / nf_modelfun_EAu (A,t,'He-4');
                
            case 'RNe'
                X(i) = nf_modelfun_EAu (A,t,'Ne-20') / nf_modelfun_EAu (A,t,'Ne-22');
                
            case 'RAr'
                X(i) = nf_modelfun_EAu (A,t,'Ar-36') / nf_modelfun_EAu (A,t,'Ar-40');
                  
            otherwise
                [c,v] = nf_atmos_gas (tracers{i},10,0,1013.25,t);                
                X(i) = v * A;

        end % switch
    end % for
end
