function X = nf_modelfun_ASW_EAu_speleowater (varargin)

% X = nf_modelfun_ASW_EAu_speleowater (T,S,P,A,M,t,tracers)
% X = nf_modelfun_ASW_EAu_speleowater (info)
%
% Returns gas amounts determined by concentrations obtained from nf_modelfun_ASW_EAu, multiplied with water amount M; see OUTPUT as described below. This is useful if water amount has significant uncertainty, which would introduce significant correlation in the errors of the concentrations calculated from gas amounts and water amount (e.g., with speleothem data). Using gas amounts and water amounts separately enables proper chi2-regression with uncorrelated errors in the input data.
%
% *** OPTION-A: X = nf_modelfun_ASW_EAu_speleowater (T,S,P,A,M,t,tracers)
%
% INPUT
% T, S, P, t: see nf_modelfun_ASW
% A: see nf_modelfun_EAu
% M: water mass (g)
%
% OUTPUT:
% X: output data (vector of entries as given in 'tracers' at input), including the water mass (M)
%
% *** OPTION-B: function X = nf_modelfun_ASW_EAu_speleowater (info)
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
            X = { 'T' 'S' 'P' 'A' 'M' 't' };

        case 'DEFAULTS'
            X    = nf_modelfun_ASW ('DEFAULTS'); % use same defaults as with nf_modelfun_ASW
            Xea  = nf_modelfun_EAu ('DEFAULTS'); % use same defaults as with nf_modelfun_ASW
            X.A  = Xea.A;
            X.M  = 1;

        case 'RANGES'
            X   = nf_modelfun_ASW ('RANGES');
            Xea = nf_modelfun_EAu ('RANGES');
            X.A = Xea.A;
            X.M.min = 0;
            X.M.max = Inf;

        case 'TRACERS'
            X = nf_modelfun_ASW_EAu ('tracers');
            X{end+1} = 'M';

        otherwise
            error (sprintf('nf_modelfun_ASW_EAu: unknown usage key ''%s''.',varargin{1}));
            
    end
else % return gas amounts
    T = varargin{1};
    S = varargin{2};
    P = varargin{3};
    A = varargin{4};
    M = varargin{5};
    t = varargin{6};
    if isempty(t)
        t = getfield (nf_modelfun_EAu('DEFAULTS'),'t');
    end
    tracers = varargin{7};
    if ( ~iscellstr(tracers) && ischar(tracers)) % convert tracers to cellstring
        tracers = cellstr (tracers);
    end
    
    % calculate sum of ASW and unfractionated excess air (amounts and isotope ratios):
    for i = 1:length(tracers)
        switch tracers{i}
                    
            case {'M'} % water amount
                X(i) = M;
                
            otherwise % gas amounts
                X(i) = nf_modelfun_ASW_EAu (T,S,P,A,t,tracers{i}) * M;
            
        end % switch
    end % for    
end
