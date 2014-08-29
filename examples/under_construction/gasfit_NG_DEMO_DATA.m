function X = nf_NG_DEMO_DATA;

% function X = nf_NG_DEMO_DATA;
%
% Returns a struct containing noble gas data of a typical groundwater sample (artificial data). Concentrations in ccSTP/g.

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