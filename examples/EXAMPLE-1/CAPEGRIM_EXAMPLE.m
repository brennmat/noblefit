% Calculate best-fit trend line describing 3He/4He ratio observed in Cape Grim Air archive tanks (1978 - 2011). See also Brennwald et al, Earth Planet. Sci. Lett., 2013, doi: 10.1016/j.epsl.2013.01.039
%
% The 3He/4He ratio is modelled to follow linear trend line: RHe(t) = RHe_0 + m*(t-2000),
% where:
%   - t: calendar year
%   - RHe: 3He/4He ratio
%   - RHe_0: 3He/4He ratio at 1.1.2000
%   - m: slope of trend line
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
% Copyright (C) 2013, 2014 Matthias S. Brennwald.
% Contact: matthias.brennwald@eawag.ch
% Further information: http://homepages.eawag.ch/~brennmat/
% *******************************************************************


% *******************************************************************
% measured data (from Brennwald et al, EPSL, 2013):
% *******************************************************************
smpl.name = 'CAPEGRIM RHe';
smpl.R1.val = 1.399E-6;	smpl.R1.err = 0.003E-6; % RHe (7.7.1978)
smpl.R2.val = 1.399E-6;	smpl.R2.err = 0.003E-6; % RHe (23.5.1984)
smpl.R3.val = 1.395E-6;	smpl.R3.err = 0.004E-6; % RHe (2.3.1993)
smpl.R4.val = 1.391E-6;	smpl.R4.err = 0.003E-6; % RHe (1.12.2004)
smpl.R5.val = 1.390E-6;	smpl.R5.err = 0.003E-6; % RHe (4.5.2011)


% *******************************************************************
% define the fit (model, fit parameters, parameter values and limits):
% *******************************************************************
mdl = 'CAPEGRIM_RHe';           	% model function
trc = {'R1','R2','R3','R4','R5'};       % tracers used in the fits: R1 (1978), R2 (1984), R3 (1993), R4 (2004), and R5 (2011)

% set which parameters are fitted and which are fixed:
fp  = [ 1 1 ]; % model parameters are RHe_0 (3He/4He on 1.1.2000) and m (slope of trend line). We want to bit both of them.

% set initial and fixed values of paramters:
p0  = [ 1.395E-6, ...        % parameter-1: RHe_0 = 1.395E-6 (3He/4He ratio)
        -1E-10 ...           % parameter-2: m = 1E-10 (slope of trend line)
];

% set minimum allowed values for the fit result:
pMin = [ 1.35E-6, ...        % parameter-1: RHe_0 >= 1.35E-6
         -1E-9 ...           % parameter-2: m >= -1E-9
];

% set maximum allowed values for the fit result:
pMax = [ 1.42E-6, ...        % parameter-1: RHe_0 <= 1.42E-6
         1E-9 ...            % parameter-2: m <= 1E-9
];

% set scaling factors (to improve numerical stability of the chi2 minimization):
pSc  = [ 1E-6,...	     % 3He/4He ratio is in the order of 1E-6
	 1E-10 ...           % slope of trend line is in the order of 1E-10
       ];


% *******************************************************************
% fit the model
% *******************************************************************
[par_val,par_err,chi2,DF,pVal] = noblefit ( mdl , smpl , trc , fp, p0 , pSc , pMin , pMax );


% *******************************************************************
% plot data and fitted trend line:
% *******************************************************************
clf;
t = [ 7/365+(7-1)/12+1978,...	% 7.7.1978
      23/365+(5-1)/12+1984,...	% 23.5.1984
      2/365+(3-1)/12+1993,...	% 2.3.1993
      1/365+(12-1)/12+2004,...	% 1.12.2004
      4/365+(5-1)/12+2011 ...	% 4.5.2011
    ];
v = [ smpl.R1.val smpl.R2.val smpl.R3.val smpl.R4.val smpl.R5.val ]; % measured values
e = [ smpl.R1.err smpl.R2.err smpl.R3.err smpl.R4.err smpl.R5.err ]; % standard errors of measured values
tt = [ 1975:10:2015 ]; mm = par_val(1) + (tt-2000)*par_val(2); % model values / trend line (tt: time, mm: 3He/4He values)
% plot measured data and fitted trend line:
h = errorbar (t,v,e,'~'); hold on % plot measured data with error bars
set (h,'marker','.','markersize',22,'linestyle','none','color','k'); % format plot
plot (tt,mm,'r'); hold off; % plot best-fit trend line
title (sprintf('3He/4He vs time -- relative change (1 sigma): (%3.2g +/- %3.2g) permil/year)',par_val(2)/par_val(1)*1000,par_err(2)/par_val(1)*1000))
