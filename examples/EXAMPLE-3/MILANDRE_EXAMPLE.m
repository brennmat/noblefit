% MILANDRE_EXAMPLE.m
%
% Script that illustrates the use of NOBLEFIT
%
% 1. Load data from file MILANDRE_SPELEO.txt (a subset of the noble gas data published in Ghadiri, E., Vogel, N., Brennwald, M. S., Maden, C., Häuselmann, A. D., Fleitmann, D., Kipfer, R. (2018). Noble gas based temperature reconstruction on a Swiss stalagmite from the last glacial–interglacial transition and its comparison with other climate records. Earth and Planetary Sciences Letters, 495, 192-201. https://doi.org/10.1016/j.epsl.2018.05.019
%
% 2. Fit a model to the data: water amount and gas amounts (air saturated water and unfractionated excess air)
%
% *******************************************************************
% This file is part of NOBLEFIT.
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
% Copyright (C) 2020 Matthias S. Brennwald.
% Contact: matthias.brennwald@eawag.ch
% Further information: http://homepages.eawag.ch/~brennmat/
% *******************************************************************


% *******************************************************************
% load data from data file
smpl = nf_read_datafile ('MILANDRE_DATA.txt');


% *******************************************************************
% define the fit (model, varialbes/parameters, limits):
mdl = 'ASW_EAu_speleowater';		% model: water amount and gas amounts (ASW and unfractionated excess air)
					% X = nf_modelfun_ASW_EAu_speleowater (T,S,P,A,M,t,tracers)
fp  = [1,0,0,1,1,0];			% fit variables: T, A and M
trc = {'Ne','Ar','Kr','Xe','M'};	% tracers used in the fits (amounts of Ne, Ar, Kr, Xe, water)

% initial / fixed parameter values
p0  = repmat ( [10,0,966,0.01,1E-4,2000] , length(smpl),1 );

% add some options for the fit:
pSc  = [];				% % scaling factors for the fitter (empty, so the fitter determines values automagically)
pMin = [0 ,0,  0,  0,1E-7,2000];	% optional: lower limits of fitted parameter values (values of fixed parameters are not used)
pMax = [40,0,Inf,0.5,1E-2,2000];	% optional: upper limits of fitted parameter values (values of fixed parameters are not used)


% *******************************************************************
% run the fit:
[par_val,par_err,chi2,DF,pVal] = noblefit ( mdl , smpl , trc , fp, p0 , pSc , pMin , pMax );


% *******************************************************************
% save results to ASCII file
nf_save_fitresults ('MILANDRE_FITRESULTS.txt',smpl,{'T (°C)','A (ccSTP/g)','M (g)'},par_val,par_err,chi2,DF,pVal,'Fit results from noblefit, fitted T-A-M fit with SPELEO model to a substet of the data in https://doi.org/10.1016/j.epsl.2018.05.019.');
