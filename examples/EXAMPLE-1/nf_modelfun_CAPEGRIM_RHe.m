function X = nf_modelfun_CAPEGRIM_RHe (varargin)

% X = nf_modelfun_CAPEGRIM_RHe (RHe_0,m)
% X = nf_modelfun_CAPEGRIM_RHe (info)
%
% Returns 3He/4He isotope ratios in air on 7.7.1978, 23.5.1984, 2.3.1993, 1.12.2004, and 4.5.2011 as calculated from a linear trend line RHe(t) = RHe_0 + (t-2000)*m,
% or usage information about the function for NOBLEGASFIT, see OUTPUT as described below.
%
% *** OPTION-A: function X = nf_modelfun_CAPEGRIM (RHe_0,m)
%
% INPUT
% RHe_0: 3He/4He ratio at 1.1.2000
% m: slope of straight line
%
% OUTPUT:
% X: output data (vector of struct with fieldnames as given in 'tracers' at input)
%
% *** OPTION-B: function X = nf_modelfun_CAPEGRIM (info)
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
            X = { 'RHe_0' , 'm' };

        case 'DEFAULTS'
            X.RHe_0 = 1.38E-6;
            X.m = 1E-10;

        case 'RANGES'
            X.RHe_0.min = 0;
            X.RHe.max = Inf;
            X.m.min = -Inf;
            X.m.max = Inf;

        case 'TRACERS'
            X = { 'R1' , 'R2' , 'R3' , 'R4' , 'R5' }; % 3He/4He ratios in 1978, 1984, 1993, 2004, and 2011

        otherwise
            error (sprintf('nf_modelfun_CAPEGRIM: unknown usage key ''%s''.',varargin{1}));

    end
else
    RHe_0 = varargin{1};
    m     = varargin{2};
    y     = varargin{3};
    if RHe_0 < 0
        warning ( 'noblefit:modelfun_params_out_of_bounds' , sprintf('nf_modelfun_CAPEGRIM: RHe_0=%g is negative. This is not physically sensible in the CAPEGRIM model...',T));
    end

    % calculate RHe values corresponding to input parameters:
    t = [ 7/365+(7-1)/12+1978,...	% 7.7.1978
	  23/365+(5-1)/12+1984,...	% 23.5.1984
          2/365+(3-1)/12+1993,...	% 2.3.1993
	  1/365+(12-1)/12+2004,...	% 1.12.2004
          4/365+(5-1)/12+2011 ...	% 4.5.2011
        ];
    X = RHe_0 + (t-2000)*m; % 3He/4He ratios at times t

end
