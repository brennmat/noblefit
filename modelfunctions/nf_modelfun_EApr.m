function X = nf_modelfun_EApr (varargin)

% X = nf_modelfun_EApr (T,A,F,t,tracers)
% X = nf_modelfun_EApr (info)
%
% Returns gas concentrations corresponding to partial re-equilibration excess air model (EAce), see OUTPUT as described below (without the ASW component). See Stute et al, Science, 269, 1995.
%
% *** OPTION-A: function X = nf_modelfun_EApr (T,A,F,t,tracers)
%
% INPUT
% T, t, tracers: see nf_modelfun_ASW
% A: amount of dry air injected per unit mass of water (ccSTP-air / g-water)
% F: fractionation parameter (scalar, dimensionless, ranging for 0 to 1). F reflects the extent of re-equilibration relative to the excess-air part of the Ne concentration (CNeEx): F = -ln (CNeEx / CNeEx(0)), where CNeEx(0) is the initial excess air concentration of Ne before re-equilibration.
%
% OUTPUT:
% X: output data (vector of struct with fieldnames as given in 'tracers' at input), see also nf_modelfun_ASW.m
%
% *** OPTION-B: function X = nf_modelfun_EApr (info)
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
            X = { 'T' , 'A' , 'F' , 't' };

        case 'DEFAULTS'
            Xasw = nf_modelfun_ASW_EAu ('DEFAULTS'); % use same defaults as with nf_modelfun_ASW_EAu
            X.T = Xasw.T;
            X.A = Xasw.A;
            X.F = 0.5;  % add fractionation parameter, F >= 0 (F=0: unfractionated EA)
            X.t = Xasw.t;

        case 'RANGES'
            Xasw = nf_modelfun_ASW_EAu ('RANGES'); % use same ranges as with nf_modelfun_ASW_EAu
            X.T = Xasw.T;
            X.A = Xasw.A;
            X.F.min = 0;    % by definition, F is not negative
            X.F.max = Inf;  % F > 0 is physically sensible
            X.t = Xasw.t;   % use same values as with nf_modelfun_ASW

        case 'TRACERS'
            X = { 'He',  'Ne' , 'Ar' , 'Kr',  'Xe' , 'N2' , 'RHe' , 'RNe' , 'RAr' , 'He_3' , 'He_4' , 'Ne_20' , 'Ne_22' , 'Ar_36' , 'Ar_40' , 'N2' , 'SF6' };

        otherwise
            error (sprintf('nf_modelfun_ASW_EApr: unknown usage key ''%s''.',varargin{1}));

    end
else

    T = varargin{1};
    A = varargin{2};
    F = varargin{3};
    t = varargin{4};
    if isempty(t)
        t = getfield (nf_modelfun_EApr('DEFAULTS'),'t');
    end
    tracers = varargin{5};
    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end
    tracers = strrep (tracers,'-','_'); % make sure entries like He-3 work out as He_3

    if A < 0
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_EApr: A=%g is negative. This is not physically sensible in the PR excess air model...',A));
    end
    
    if F < 0
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_EApr: F=%g is negative. This is not physically sensible in the PR model...',F));
    end
    
    X = repmat (NA,1,length(tracers));   
    
    for i = 1:length(tracers) % calculate according to Eqn (10) in Aeschbach-Hertig et al., WRR, 2002
        switch tracers{i}
        
            case 'RHe'
                X(i) = nf_modelfun_EApr (T,A,F,t,'He-3') / nf_modelfun_EApr (T,A,F,t,'He-4');

            case 'RNe'
                X(i) = nf_modelfun_EApr (T,A,F,t,'Ne-20') / nf_modelfun_EApr (T,A,F,t,'Ne-22');

            case 'RAr'
                X(i) = nf_modelfun_EApr (T,A,F,t,'Ar-36') / nf_modelfun_EApr (T,A,F,t,'Ar-40');
                  
            otherwise
                X(i) = nf_modelfun_EAu (A,t,tracers{i}); % get unfractionated excess air
                % calculate fractionation:
                if ( A > 0 && F > 0 ) % only worry about fractionation if A and F are both non-zero
                    if ~strcmp (upper(tracers{i}),'NE')
                        [u1,u2,Di] = nf_atmos_gas (tracers{i},T,0,1013.25,t);    % Di: diffusion coefficient of tracers{i}
                        [u1,u2,DNe] = nf_atmos_gas ('Ne',T,0,1013.25,t);         % DNe: Ne diffusion coefficient
                        X(i) = X(i) * exp (-F * Di/DNe); % add fractionation part due to partial re-equilibration
                    else % don't waste time with calculating diffusion coefficients and things:
                        X(i) = X(i) * exp (-F); % add fractionation part due to partial re-equilibration
                    end
                end % if A > 0 and F > 0
        end % switch
    end % for
end
