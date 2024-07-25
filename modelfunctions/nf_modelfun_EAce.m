function X = nf_modelfun_EAce (varargin)

% X = nf_modelfun_EAce (T,S,P,A,F,t,tracers)
% X = nf_modelfun_EAce (info)
%
% Returns gas concentrations corresponding to closed-system equilibrium excess air (EAce), see OUTPUT as described below (without the ASW component). See Aeschbach-Hertig et al., Nature, 405, 2000.
%
% *** OPTION-A: function X = nf_modelfun_EAce (T,S,P,A,F,tracers)
%
% INPUT
% T, S, P, t, tracers: see nf_modelfun_ASW
% A: initial amount of entrapped air per unit mass of water (ccSTP-air / g-water)
% F: fractionation parameter (scalar, dimensionless, ranging for 0 to 1)
%
% OUTPUT:
% X: output data (vector of entries as given in 'tracers' at input), see also nf_modelfun_ASW.m
%
% *** OPTION-B: function X = nf_modelfun_EAce (info)
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
%   info = 'DEFAULTS' returns struct with default values that may be used as initial values for fitting.
%   info = 'RANGES' returns struct with fields min / max values
%   info = 'TRACERS' returns cell string with tracer names
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
            X = { 'T' 'S' 'P' 'A' 'F' };

        case 'DEFAULTS'
            X    = nf_modelfun_ASW_EAu ('DEFAULTS'); % use same defaults as with nf_modelfun_ASW_EAu
            X.F  = 0.2;  % add fractionation parameter, 0 <= F <= 1 (F=0: unfractionated EA)

        case 'RANGES'
            X    = nf_modelfun_ASW_EAu ('RANGES'); % use same ranges as with nf_modelfun_ASW_EAu
            X.F.min = 0; % by definition, F is not negative
            X.F.max = 1; % by definition, F > corresponds to a degassing model (which is not the idea here)

        case 'TRACERS'
            X = { 'He',  'Ne' , 'Ar' , 'Kr',  'Xe' , 'N2' , 'RHe' , 'RNe' , 'RAr' , 'He_3' , 'He_4' , 'Ne_20' , 'Ne_22' , 'Ar_36' , 'Ar_40' , 'Kr_84' , 'Kr_86' , 'Xe_136' , 'N2' , 'SF6' };

        otherwise
            error (sprintf('nf_modelfun_ASW_EAce: unknown usage key ''%s''.',varargin{1}));

    end
else

    T = varargin{1};
    S = varargin{2};
    P = varargin{3};
    A = varargin{4};
    F = varargin{5};
    t = varargin{6};
    if isempty(t)
        t = getfield (nf_modelfun_Ece('DEFAULTS'),'t');
    end
    tracers = varargin{7};
    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end
    tracers = strrep (tracers,'-','_'); % make sure entries like He-3 work out as He_3


    if A < 0
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_EAce: A=%g is negative. This is not physically sensible in the CE excess air model...',A));
    end
    
    if F > 1
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_EAce: F=%g is larger than 1. This is not physically sensible in the CE model...',F) );
    elseif F < 0
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_EAce: F=%g is negative. This is not physically sensible in the CE model...',F));
    end
    
    X = repmat (NA,1,length(tracers));   
    
    for i = 1:length(tracers) % calculate according to Eqn (3) in Aeschbach-Hertig et al., Nature, 405, 2000
        switch tracers{i}
        
            case 'RHe'
                X(i) = nf_modelfun_EAce (T,S,P,A,F,t,'He-3') / nf_modelfun_EAce (T,S,P,A,F,t,'He-4');

            case 'RNe'
                X(i) = nf_modelfun_EAce (T,S,P,A,F,t,'Ne-20') / nf_modelfun_EAce (T,S,P,A,F,t,'Ne-22');

            case 'RAr'
                X(i) = nf_modelfun_EAce (T,S,P,A,F,t,'Ar-36') / nf_modelfun_EAce (T,S,P,A,F,t,'Ar-40');
                  
            otherwise
                [C,z] = nf_atmos_gas (tracers{i},T,S,P,t);                
                X(i) = (1-F)*A*z / (1+F*A*z/C) ; 
                                
        end % switch
    end % for
end
