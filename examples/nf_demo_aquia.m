% This script (nf_demo_aquia) illustrates the use of NOBLEFIT to reproduce the fit of noble gas data obtained in the Aquia Aquifer (see also W. Aeschbach-Hertig, M. Stute, J. Clark, R. Reuter, and P. Schlosser. A paleotemperature record derived from dissolved noble gases in groundwater of the Aquia Aquifer (Maryland, USA). Geochim. Cosmochim. Acta, 66(5):797–817, 2002.):
%
% 1. Load data from file AQUIA_GW.txt (a subset of the noble gas data published in W. Aeschbach-Hertig, M. Stute, J. Clark, R. Reuter, and P. Schlosser. A paleotemperature record derived from dissolved noble gases in groundwater of the Aquia Aquifer (Maryland, USA). Geochim. Cosmochim. Acta, 66(5):797–817, 2002.)
%
% 2. Fit a model to the data: ASW(T) and EA(T,A,F) according to the CE model.
%
%
% *******************************************************************
% This file is part of NOBLEFIT.
% 
% NOBLEFIT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% NOBLEFIT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with NOBLEFIT.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) 2013 Matthias S. Brennwald.
% Contact: matthias.brennwald@eawag.ch
% Further information: http://homepages.eawag.ch/~brennmat/
% *******************************************************************


% *******************************************************************
% load measured noble gas data from file:
% *******************************************************************
options.replace_zeros           = NA; % replace zero values in data file by NA
options.filter_InvIsotopeRatios = 1;  % try to replace isotope ratios in the form heavy/light by the inverse (light/heavy)
smpl = nf_read_datafile ('AQUIA_GW.txt',options);


% *******************************************************************
% define the fit (model, fit parameters, parameter values and limits):
% *******************************************************************
mdl = 'ASW_EAce';               % model: ASW and closed-system excess air
trc = {'Ne','Ar','Kr','Xe'};    % tracers used in the fits: Ne, Ar, Kr, and Xe

% set which parameters are fitted and which are fixed:
fp  = [ 1, ...  % parameter-1: T, is fitted
        0, ...  % parameter-2: S, fixed
        0, ...  % parameter-3: p, fixed
        1, ...  % parameter-4: A, fitted
        1, ...  % parameter-5: F, fitted
        0 ...   % parameter-6: t, fixed
];

% set initial and fixed values of paramters (see also Aeschbach et al., GCA, 2002):
p0  = [ 8, ...              % parameter-1: T = 8 deg.C (initial value)
        0, ...              % parameter-2: S = 0 permil (fixed value)
        0.944*1013.25, ...  % parameter-3: p = 0.994 atm = 0.944*1013.25 hPa (fixed value)
        1E-3, ...           % parameter-4: A = 1E-3 ccSTP/g (initial value)
        0, ...              % parameter-5: F = 0 (initial value)
        2000 ...            % parameter-6: t = 2000 A.D. (fixed value)
];

% set minimum allowed values:
pMin = [ 0, ...              % parameter-1: T >= 0 (fitted value)
         0, ...              % parameter-2: S >= 0 permil (fixed value)
         0.944*1013.25, ...  % parameter-3: p >= 0.994 atm = 0.944*1013.25 hPa (fixed value)
         0, ...              % parameter-4: A >= 1E-3 ccSTP/g (fitted value)
         0, ...              % parameter-5: F >= 0 (fitted value)
         2000 ...            % parameter-6: t = 2000 A.D. (fixed value)
];

% set maximum allowed values:
pMax = [ 24, ...             % parameter-1: T <= 25 (fitted value)
         0, ...              % parameter-2: S <= 0 permil (fixed value)
         0.944*1013.25, ...  % parameter-3: p <= 0.994 atm = 0.944*1013.25 hPa (fixed value)
         0.1, ...              % parameter-4: A <= 0.1 ccSTP/g (fitted value)
         1, ...              % parameter-5: F <= 1 (fitted value)
         2000 ...            % parameter-6: t = 2000 A.D. (fixed value)
];

% set scaling factors
pSc  = []; % don't specify scaling factors (the fitter determines values automagically)


% *******************************************************************
% fit the model
% *******************************************************************
[par_val,par_err,chi2,DF,pVal] = noblefit ( mdl , smpl , trc , fp, p0 , pSc , pMin , pMax );
