function G = nf_objfun (PF,P0,PFmin,PFmax,mdl,X_val,X_err)

% G = nf_objfun (PF,P0,PFmin,PFmax,mdl,X_val,X_err)
%
% Objective function for mimization in parameter regression. Determines the modelled values using the model for the given parameter values, and calculates the chi^2 value from the difference to the data values. If called with fit parameter values that exceed the limits in PFmin and PFmax, the chi^2 value will be calculated such that the chi^2 minimizer will look elsewhere for a 'better' optimum (see commented code for details). Note that many techniques for fitting with constrained parameter ranges rely on mapping the "infinite chi2 surface" to the allowed parameter range. This distorts the chi2 surface, and while the minimum chi2 value will correspond to the correct best-fit paramter values, the "distorted" chi2 value cannot be statistically interpreted in the same way as the "undistorted" chi2. nf_objfun.m therefore avoids distorting the chi2 surface in the allowed parameter range, and the resulting chi2 value can be interpreted using the conventionel chi2 statistics.
%
% INPUT:
% PF: vector of fitted model parameter values
% P0: vector of constant model parameter values
% PFmin: vector of minimum values allowed for the fitted parameters
% PFmax: vector of maximum values allowed for the fitted parameters
% mdl: string containing call to the model function using PF and P0
% X_val: values of observed/measured data (rows correspond to samples, columns correspond to one tracers)
% X_err: standard errors of X_val
%
% OUTPUT:
% G: chi^2 value
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
% Copyright (C) 2014 Matthias S. Brennwald, matthias.brennwald@eawag.ch
% Additional information about NOBLEFIT is available at http://homepages.eawag.ch/~brennmat/
% Please contribute if you find this software useful.
% *******************************************************************

% global last_chi2; % used to deal with parameter values outside "the allowed range"

% THIS IS DONE IN THE DEFINITION OF THE INLINE FUNCTION CALL TO nf_objfun --- PF = (PF+PF0) .* PFnorm; % undo parameter normalisaiton before calling model function

fflush (stdout); % flush out any messages produced by the model function



% Check if PF values are within allowed limits:

if any ( k = (PF < PFmin) ) % some parameter values are lower than allowed.
    W1 = warning ('query','nf_atmos_gas_timerange'); warning ('off','nf_atmos_gas_timerange'); 
    W2 = warning ('query','nf_atmos_gas_henry_zerodivision'); warning ('off','nf_atmos_gas_henry_zerodivision'); 
    M_val = eval (mdl); % avaluate model in the forbidden area
    warning (W1.state,'nf_atmos_gas_timerange');
    warning (W2.state,'nf_atmos_gas_henry_zerodivision');
    Gx = sum ( ((X_val-M_val) ./ X_err).^2 ); % chi2 in the forbidden area
    dP = PFmin(k) - PF(k); % the delta relative to the limit values
    l0   = sqrt(sum(dP.^2)); % length of dP (remember for later)
    % jump aver the limit value(s) from the forbidden area to the allowed area:
    while any ( kk = (PFmin(k)+dP > PFmax(k)) ) % we'd jump too far...
        dP(kk) = dP(kk) / 2; % make the jump smaller (until we end up in the allowed area)
    end
    PF(k) = PFmin(k) + dP; % jump into the allowed area
    M_val = eval (mdl); % evaluate model in the allowed area
    Ga = sum ( ((X_val-M_val) ./ X_err).^2 ); % chi2 in the allowed area
    PF(k) = PFmin(k); % this on the border to the allowed area
    M_val = eval (mdl); % avaluate model on the border
    Gb = sum ( ((X_val-M_val) ./ X_err).^2 ); % chi2 on the border
    Gf = Gb + abs(Ga-Gb) * l0/sqrt(sum(dP.^2)); % fake the value in the forbidden area by mirroring the chi2 curve at the limit value (by linear extrapolation)

    G = max(Gx,Gf); % take the larger of the two
    
elseif any ( k = (PF > PFmax) ) % some parameter values are larger than allowed.
    W1 = warning ('query','nf_atmos_gas_timerange'); warning ('off','nf_atmos_gas_timerange'); 
    W2 = warning ('query','nf_atmos_gas_henry_zerodivision'); warning ('off','nf_atmos_gas_henry_zerodivision'); 
    M_val = eval (mdl); % avaluate model in the forbidden area
    warning (W1.state,'nf_atmos_gas_timerange');
    warning (W2.state,'nf_atmos_gas_henry_zerodivision');
    Gx = sum ( ((X_val-M_val) ./ X_err).^2 ); % chi2 in the forbidden area

    dP = PF(k) - PFmax(k); % the delta relative to the limit values
    l0   = sqrt(sum(dP.^2)); % length of dP (remember for later)
    % jump aver the limit value(s) from the forbidden area to the allowed area:
    while any ( kk = (PFmax(k)-dP < PFmin(k)) ) % we'd jump too far...
        dP(kk) = dP(kk) / 2; % make the jump smaller (until we end up in the allowed area)
    end
    PF(k) = PFmax(k) - dP; % jump into the allowed area
    M_val = eval (mdl); % evaluate model in the allowed area
    Ga = sum ( ((X_val-M_val) ./ X_err).^2 ); % chi2 in the allowed area
    PF(k) = PFmin(k); % this on the border to the allowed area
    M_val = eval (mdl); % avaluate model on the border
    Gb = sum ( ((X_val-M_val) ./ X_err).^2 ); % chi2 on the border
    Gf = Gb + abs(Ga-Gb) * l0/sqrt(sum(dP.^2)); % fake the value in the forbidden area by mirroring the chi2 curve at the limit value (by linear extrapolation)

    G = max(Gx,Gf); % take the larger of the two
    
else % Parameter value(s) in the "allowed range". Go ahead and calculate chi2:
    M_val = eval (mdl);
    G = sum ( ((X_val-M_val) ./ X_err).^2 );
    
end
