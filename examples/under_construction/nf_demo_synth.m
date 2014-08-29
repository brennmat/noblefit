ID = round (1000*rand(1))

% diary (tilde_expand(sprintf('~/ngf_demo_synth_%i.log',ID)));

% tracers included in the synthetic data set:
trc = { 'He' 'Ne','Ar','Kr','Xe' };    % tracers included in the synthetic data

% parameters used to calculate synthetic data:

N = 100; % number of synthetic 'samples'


parname = { 'T' 'S' 'p' 'A' 'F' };

iFit =  [ 1 0 0 1 1 ]; % which parameters should be varied and and then fitted (T S P A F) ?
T_min = 0;
T_max = 30;
S_min = 0;
S_max = 0;
p_min = 1013.25;
p_max = 1013.25;
A_min = 1E-4;
A_max = 1E-2;
F_min = 1E-2;
F_max = 1;

for i = 1:length(iFit)
    if iFit(i) % variable parameter
        eval (sprintf('%s = %s_min + rand(1,N)*(%s_max-%s_min);',parname{i},parname{i},parname{i},parname{i}));
    else
        eval (sprintf('%s = repmat((%s_min + %s_max)/2,1,N);',parname{i},parname{i},parname{i}));
    end
end

% calculate synthetic data (incl. random errors with Gaussion distribution)
disp ('Calculating synthetic data...'); fflush (stdout);
rel_err = 0.01; % relative standard error of concentrations
k = 1;
conc = [];
for i = 1:N
        
    u_val = ngf_modelfun_ASW_EApr ( T(i) , S(i) , p(i) , A(i) , F(i) , trc ); % concentration values   
    u_err = rel_err * u_val; % standard error

    u_val = u_val + u_err .* randn(size(u_err)); % add random errors to calculated concentration value (with Gaussian distribution)
    
    for i_trc = 1:length(trc)
        eval (sprintf ('u.%s.val = u_val(i_trc);',trc{i_trc}));
        eval (sprintf ('u.%s.err = u_err(i_trc);',trc{i_trc}));
    end
    eval ( 'u.T = T(i);' );
    eval ( 'u.S = S(i);' );
    eval ( 'u.p = p(i);' );
    eval ( 'u.A = A(i);' );
    eval ( 'u.F = F(i);' );

    conc = [ conc ; u ];
        
    k = k+1;
end

% *******************************************************************
% define the fit (model, parameters, limits):
mdl  = 'ASW_EApr';                   % model: ASW and excess air (partial re-equilibration model)
p0   = [mean(T),mean(S),mean(p),mean(A),mean(F)];   % initial / fixed values for the fits
pSc  = [5 5 1000 1E-3 0.5 ];         % scaling factors for the fitter (empty, so the fitter determines values automagically)
pMin = [0,0,0,0,0];                  % lower limits of fitted parameter values (values of fixed parameters are not used)
pMax = [];                           % upper limits of fitted parameter values (values of fixed parameters are not used) 


% *******************************************************************
        
[val,err,chi2,DF,pVal] = noblegasfit ( mdl , conc , trc , iFit, p0 , pSc , pMin , pMax ); % fit parameters that are variable in synthetic data

% save results do disk
f = tilde_expand (sprintf('~/ngf_synth_fitresults_%i.mat',ID));
save ( f, 'conc' , 'mdl' , 'val' , 'err' , 'chi2' , 'DF' , 'pVal' , 'iFit' , 'p0' , 'pSc' , 'pMin' , 'pMax' , 'parname' );

diary off