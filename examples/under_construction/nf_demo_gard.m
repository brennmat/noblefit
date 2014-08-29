% ngf_demo_gard
%
% Script that illustrates the use of NOBLEGASFIT.
%
% 1. Load data from file GARD_GW.txt (noble gas data from Gardermoen aquifer, Norway)
%
% 2. Fit a model to the data: ASW(T) and unfractionated excess air EA(T,A).
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
% Copyright (C) 2013 Matthias S. Brennwald.
% Contact: matthias.brennwald@eawag.ch
% Further information: http://homepages.eawag.ch/~brennmat/
% *******************************************************************


% *******************************************************************
% load example data from data file
options.replace_zeros           = NA; % replace zero values in data file by NA
options.filter_InvIsotopeRatios = 1;  % try to replace isotope ratios in the form heavy/light by the inverse (light/heavy)
smpl = ngf_read_datafile ('GARD_GW.txt',options);


% *******************************************************************
% define the fit (model, parameters):
mdl = 'ASW_EAu';                % model: ASW and unfractionated excess air
trc = {'Ne','Ar','Kr','Xe'};    % tracers used in the fits: Ne, Ar, Kr, and Xe
fp  = [1,0,0,1,0];                % fit first and fourth parameters of model function (temperature T and excess air concentration A)
p0  = [5,0,988.53,1E-3,2000];        % initial / fixed values for the fits (using the same values for all samples). T-initial = 5 deg.C, S = 0 g/kg (fixed), p = 988.53 hPa (fixed), A-initial = 1E-3 ccSTP/g, year = 2000 (irrelevant for non-transient atmospheric gases)


% *******************************************************************
% fit all samples in smpl (parameter for fitting is done automagically, fit parameters have no upper or lower bounds):
[par_val,par_err,chi2,DF,pVal] = noblegasfit ( mdl , smpl , trc , fp, p0 );


% *******************************************************************
% add some options for the fit:
pSc  = [1,1,1E3,1E-2,1];          % optional: scaling factors for the fitter (so that parameter values are scaled to similar numerical ranges during fitting)
pMin = [0,0,0,0,2000];               % optional: lower limits of fitted parameter values (values of fixed parameters are not used)
pMax = [25,0,Inf,0.1,2000];          % optional: upper limits of fitted parameter values (values of fixed parameters are not used) 


% *******************************************************************
% repeat the fitting, this time with the additional options:
[par_val,par_err,chi2,DF,pVal] = noblegasfit ( mdl , smpl , trc , fp, p0 , pSc , pMin , pMax );
