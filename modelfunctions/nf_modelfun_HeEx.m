function X = nf_modelfun_HeEx (varargin)

% X = nf_modelfun_HeEx (H,R,TR,tracers)
% X = nf_modelfun_HeEx (info)
%
% Returns concentrations gas resulting from non-atmospheric He isotope sources (terrigenic He and tritiogenic He-3), see OUTPUT as described below.
%
% *** OPTION-A: function X = nf_modelfun_HeEx (H,R,TR,tracers)
%
% INPUT
% H: amount of terrigenic He accumulated per unit mass of water (H = He-3 + He-4), (ccSTP-He / g-water)
% R: He-3 / He-4 isotope ratio of terrigenic He
% TR: amount of tritiogenic He-3 accumulated per unit mass of water, (ccSTP-He-3 / g-water)
% tracers: see nf_modelfun_ASW
%
% OUTPUT:
% X: output data (vector of struct with fieldnames as given in 'tracers' at input), see also nf_modelfun_ASW.m
%
% *** OPTION-B: function X = nf_modelfun_HeEx (info)
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
% Copyright (C) 2019 Matthias S. Brennwald, matthias.brennwald@eawag.ch
% Additional information about NOBLEGASFIT is available at http://homepages.eawag.ch/~brennmat/
% Please contribute if you find this software useful.
% *******************************************************************

if nargin == 1
    switch upper(varargin{1})
    
        case 'PARAMS'
            X = { 'H' , 'R' , 'TR' };

        case 'DEFAULTS'
            X.H = 1E-10; % 1E-10 ccSTP-He-rad / g-water
            X.R = 2E-8; % R = 2E-8, typical radiogenic 3He/4He ratio
            X.TR = 1E-15; % TR = 1E-15 ccSTP-He3-trit / g-water
            
        case 'RANGES'
            X.H.min = 0;
            X.H.max = Inf;
            X.R.min = 0;
            X.R.max = Inf;
            X.TR.min = 0;
            X.TR.max = Inf;

        case 'TRACERS'
            X = { 'He',  'Ne' , 'Ar' , 'Kr',  'Xe' , 'N2' , 'RHe' , 'RNe' , 'RAr' , 'He_3' , 'He_4' , 'Ne_20' , 'Ne_22' , 'Ar_36' , 'Ar_40' , 'Kr_84' , 'Kr_86' , 'Xe_136' , 'N2' , 'SF6' };

        otherwise
            error (sprintf('nf_modelfun_HeEx: unknown usage key ''%s''.',varargin{1}));

    end
else
    H  = varargin{1};
    R  = varargin{2};
    TR = varargin{3};
    tracers = varargin{4};
    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end
    tracers = strrep (tracers,'-','_'); % make sure entries like He-3 work out as He_3

    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end

    if H < 0
        warning ('noblefit:modelfun_params_out_of_bounds',sprintf('nf_modelfun_HeEx: H=%g is negative. This is not physically sensible...',H));
    end
    if R < 0
        warning ('noblefit:modelfun_params_out_of_bounds',sprintf('nf_modelfun_HeEx: R=%g is negative. This is not physically sensible...',R));
    end
    if TR < 0
        warning ('noblefit:modelfun_params_out_of_bounds',sprintf('nf_modelfun_HeEx: TR=%g is negative. This is not physically sensible...',TR));
    end

    X = repmat (NA,1,length(tracers));

    for i = 1:length(tracers)
        switch tracers{i}
        
            case 'RHe'
                X(i) = R;
                
            case 'He'
                X(i) = H + TR; % sum of radiogenic He and tritiogenic He

            case 'He_4'
                X(i) = H / (1+R) ; % radiogenic He-4
                
            case 'He_3'
                X(i) = R * ( H / (1+R) ) + TR; % radiogenic He-3 + tritigenic He-3
                  
            otherwise                
                X(i) = 0; % all non-Helium components are zero

        end % switch
    end % for
end
