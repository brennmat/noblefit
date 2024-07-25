function X = nf_modelfun_ASW_EAu (varargin)

% X = nf_modelfun_ASW_EAu (T,S,P,A,t,tracers)
% X = nf_modelfun_ASW_EAu (info)
%
% Returns gas concentrations given by sum of air saturated water (ASW) and unfractionation excess air (EAu), see OUTPUT as described below.
%
% *** OPTION-A: X = nf_modelfun_ASW_EAu (T,S,P,A,t,tracers)
%
% INPUT
% T, S, P, t, tracers: see nf_modelfun_ASW
% A: see nf_modelfun_EAu
%
% OUTPUT:
% X: output data (vector of entries as given in 'tracers' at input), see also nf_modelfun_ASW.m
%
% *** OPTION-B: function X = nf_modelfun_ASW_EAu (info)
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

if nargin == 1 % return info
    switch upper(varargin{1})

        case 'PARAMS'
            X = { 'T' 'S' 'P' 'A' 't' };

        case 'DEFAULTS'
            X    = nf_modelfun_ASW ('DEFAULTS'); % use same defaults as with nf_modelfun_ASW
            Xea  = nf_modelfun_EAu ('DEFAULTS'); % use same defaults as with nf_modelfun_ASW
            X.A  = Xea.A;

        case 'RANGES'
            X   = nf_modelfun_ASW ('RANGES');
            Xea = nf_modelfun_EAu ('RANGES');
            X.A = Xea.A;

        case 'TRACERS'
            tASW = nf_modelfun_ASW ('TRACERS');
            tEA  = nf_modelfun_EAu ('TRACERS');
            X    = intersect (tASW,tEA); % tracers that are included in both model parts

        otherwise
            error (sprintf('nf_modelfun_ASW_EAu: unknown usage key ''%s''.',varargin{1}));
            
    end
else % return concentrations
    T = varargin{1};
    S = varargin{2};
    P = varargin{3};
    A = varargin{4};
    t = varargin{5};
    if isempty(t)
        t = getfield (nf_modelfun_EAu('DEFAULTS'),'t');
    end
    tracers = varargin{6};
    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end
    tracers = strrep (tracers,'-','_'); % make sure entries like He-3 work out as He_3
    
    % calculate sum of ASW and unfractionated excess air (concentrations and isotope ratios):
    for i = 1:length(tracers)
        switch tracers{i}
        
            case 'RHe'            
                X(i) = nf_modelfun_ASW_EAu (T,S,P,A,t,'He-3') / nf_modelfun_ASW_EAu (T,S,P,A,t,'He-4');

            case 'RNe'
                X(i) = nf_modelfun_ASW_EAu (T,S,P,A,t,'Ne-20') / nf_modelfun_ASW_EAu (T,S,P,A,t,'Ne-22');
            
            case 'RAr'
                X(i) = nf_modelfun_ASW_EAu (T,S,P,A,t,'Ar-36') / nf_modelfun_ASW_EAu (T,S,P,A,t,'Ar-40');
            
            otherwise % concentrations
                X(i)  = nf_modelfun_ASW (T,S,P,t,tracers{i}) + nf_modelfun_EAu (A,t,tracers{i});                
                
        end % switch
    end % for    
end
