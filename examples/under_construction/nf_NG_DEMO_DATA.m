function X = ngf_NG_DEMO_DATA;

% function X = ngf_NG_DEMO_DATA;
%
% Returns a struct containing noble gas data of a typical groundwater sample (artificial data). Concentrations in ccSTP/g.
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

T = 7;   % water temperature in deg.C
A = 1E-3;       % excess air concentration in ccSTP-air/g
S = 0;
P = exp (-500/8300) * 1013.25;

[X.He.val,v] = atmos_gas_ccSTP ('He',T,S,P);
X.He.val = X.He.val + v*A;
X.He.err = 0.015 * X.He.val;

[X.Ne.val,v] = atmos_gas_ccSTP ('Ne',T,S,P);
X.Ne.val = X.Ne.val + v*A;
X.Ne.err = 0.015 * X.Ne.val;

[X.Ar.val,v] = atmos_gas_ccSTP ('Ar',T,S,P);
X.Ar.val = X.Ar.val + v*A;
X.Ar.err = 0.015 * X.Ar.val;

[X.Kr.val,v] = atmos_gas_ccSTP ('Kr',T,S,P);
X.Kr.val = X.Kr.val + v*A;
X.Kr.err = 0.015 * X.Kr.val;

[X.Xe.val,v] = atmos_gas_ccSTP ('Xe',T,S,P);
X.Xe.val = X.Xe.val + v*A;
X.Xe.err = 0.015 * X.Xe.val;