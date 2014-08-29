function [par_val,par_err,chi2,DF,pVal,cov,res] = noblefit (model,tracer_data,tracers,par_usage,par0,par_norm,par_min,par_max);

% [par_val,par_err,chi2,DF,pVal,cov,res] = noblefit (model,tracer_data,tracers,par_usage,par0,par_norm,par_range,par_min,par_max);
%
% Frontend to fit a given model to observed data using chi^2 regression.
%
% INPUT:
% model: model name (string). More specifically, this is the name of the function that returnscan either be the name of one of the standard models provided with gasfit, or a custom model function provided by the user (function name must be on the search path).
% tracer_data: either the name of an ASCII file containing the data (with column headers that correspond to the tracer names) or a struct variable containing the observed data (either a single struct corresponding to one sample, or a vector of structs for more than one sample). The tracer names must correspond to those in the model function.
% tracers: a cell string containing the names of the tracers that are to be used in the fit (names must correspond to those used by the model function).
% par_usage: vector indicating usage of the model parameter values. par_usage is of the same format as the vector used as the input argument of the model function. Values are as follows:
%   - par_usage(i) = 0: The i-th model parameter is used as a constant in the minimization problem, i.e., it is not optimized during fitting of the model to the data.
%   - par_usage(i) = 1: the i-th parameter is optimized to obtain the best model fit for each individual sample.
%   - other values may be implemented in a later version (e.g., for ensemble fits of a model parameter to an ensemble of data from multiple samples)
% par0: parameter values used for the regression (vector of the same format as used for the input argument of the model function). Depending on the parameter usage (see par_usage), the values are used as fixed values or initial values used in the minimization problem.
% par_norm (optional): typical scale of variation of the parameter values, used to normalise fitting parameters during model fitting (vector or matrix of sime size as par0, values used for fit parameters must not be zero). Note that the scaling factors reflect the RANGE OF VARIATION, not the absolute value. For instance, if infiltration date is somewhere between 1950 and 2013, a suitable scaling factor would be 10, not 1000.
% par_min (optional): min allowed parameter values for fit (vector or matrix of same size as par0, only values for fitted parameters will be used). Values may be -Inf to indicate no limits. Note: fit results may sightly exceed the limits, if the best fit value is close to the limit. This is due to the way the limits are treated in the fitting routine. The effect should be small, but please make sure the values are ok for you.
% par_max (optional): max allowed parameter values for fit (vector or matrix of same size as par0, only values for fitted parameters will be used). Values may be Inf to indicate no limits. Note: fit results may sightly exceed the limits, if the best fit value is close to the limit. This is due to the way the limits are treated in the fitting routine. The effect should be small, but please make sure the values are ok for you.
%
% OUTPUT:
% par_val: best-fit estimates of the parameter values (vector of the same format as used for the input argument of the model function).
% par_err: standard errors of par_val (vector).
% chi2: chi2 of fit
% DF: degrees of freedom of fit ( = number of data points - number of fitted model parameters)
% pVal: p-value of chi2, as given from cumulative chi2-density function with DF degrees of freedom: pVal = 1 - chi2cdf(chi2,DF)
% cov: covariance matrices of the fits (cell array of matrices)
% res: error-weighted residuals of observed values (X) of sample k relative to modeled values (M), normalised by the standard error (E) of the observed values, i.e.: res(k,i) = (X_i - M_i) / E_i, where i is an index to the corresponding tracer (res is a matrix, each row corresponds to one sample, columns correspond to 'tracers')
%
% EXAMPLES:
% see files in Examples folder.
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


% ****** CHECK USER INPUT AND PRINT SUMMARY ******

disp (''); disp (''); 
disp ('***************************************************************')
disp ('*** NOBLEFIT: CHECKING INPUT AND PREPARING MODEL FITTING ***')
disp ('***************************************************************')

model_fun       = sprintf('nf_modelfun_%s',model); % name of the model function
model_fun_file  = which (model_fun); % function file of the model

if isempty(model_fun_file)
    error (sprintf('noblefit: unknown model function %s.',model_fun))
end

if ~iscellstr(tracers)
    if ischar(tracers)
        tracers = cellstr (tracers);
    end
end

model_tracers  = eval (sprintf('%s(''TRACERS'');',model_fun)); % list of tracers included in the model
model_pars     = eval (sprintf('%s(''PARAMS'');',model_fun)); % list of input parameters of model function (cell string)
Nsmpl          = length (tracer_data);  % number of samples
Ntrcr          = length (tracers); % number of tracers used in the chi2 minimisation
NFpar          = sum (par_usage == 1);  % number of fit parameters
cov            = cell (Nsmpl,1); % cell array of covariance matrices
res            = repmat (NA,Nsmpl,Ntrcr); % residuals matrix


if isempty (model_fun_file) % produce an error (and abort)
    error (sprintf('gasfit.m error: model function is not defined for model ''%s'' (expected function name: %s).',model,model_fun))
end

% check consistency of input

% check tracers:
if Ntrcr < 1
    error ('noblefit: ''tracers'' is empty. Cannot fit anything without including any tracers.');
end
for i = 1:Ntrcr % construct string containing tracers included in the fit, and check agreement of tracer names
    if ~any(strcmp(tracers{i},model_tracers))
        error (sprintf('noblefit: the model function ''%s'' does not provide a tracer named ''%s''. Please check the names of the tracers included in the model function.',model_fun,tracers{i}));
    end
    if i == 1
        trc = tracers{i};
    else
        trc = sprintf ('%s, %s',trc,tracers{i});
    end
end

% check par0
if size(par0,2) ~= length(model_pars)
    error (sprintf('noblefit: number of columns in ''par0'' (%i) must agree with the number of model parameters (%i) of the model function (%s).',size(par0,2),length(model_pars),model_fun))
end
if size (par0,1) ~= Nsmpl
    if size(par0,1) == 1
        warning (sprintf('par0 contains only one row, but there are %i samples. Assuming same values of par0 for all samples...',Nsmpl));
        par0 = repmat (par0,Nsmpl,1);
    else
        error (sprintf('Number of rows in par0 (%i) does not match number of samples (%i).',size(par0,1),Nsmpl));
    end
end

% check par_usage
if size(par_usage,2) ~= length(model_pars)
    error (sprintf('noblefit: number of columns in ''par_usage'' (%i) must agree with the number of model parameters (%i) of the model function (%s).',size(par_usage,2),length(model_pars),model_fun))
end

% check par_norm
if exist ('par_norm','var')
    if isempty (par_norm) % an empty vector was specified
        par_norm = repmat (NA,size(par0));
    elseif size(par_norm,2) ~= length(model_pars)
        error (sprintf('noblefit: number of columns in ''par_norm'' (%i) must agree with the number of model parameters (%i) of the model function (%s).',size(par_norm,2),length(model_pars),model_fun))
    end
    if size (par_norm,1) ~= Nsmpl
        if size(par_norm,1) == 1
            warning (sprintf('par_orm contains only one row, but there are %i samples. Assuming same values of par_norm for all samples...',Nsmpl));
            par_norm = repmat (par_norm,Nsmpl,1);
        else
            error (sprintf('Number of rows in par0 (%i) does not match number of samples (%i).',size(par0,1),Nsmpl));
        end
    end
else % try to guess suitable values for par_norm
    disp ('Parameter normalisation for fitting not specified! Will try to determine suitable values later (see below)...')
    par_norm = repmat (NA,size(par0));
end

% check par_min
if ~exist ('par_min','var')
    par_min = repmat (-Inf,size(par0));
end
if isempty (par_min) % an empty vector was specified
    par_min = repmat (-Inf,size(par0));
elseif size(par_min,2) ~= length(model_pars)
    error (sprintf('noblefit: number of columns in ''par_min'' (%i) must agree with the number of model parameters (%i) of the model function (%s).',size(par_min,2),length(model_pars),model_fun))
end
if size (par_min,1) ~= Nsmpl
    if size(par_min,1) == 1
        warning (sprintf('par_min contains only one row, but there are %i samples. Assuming same values of par_min for all samples...',Nsmpl));
        par_min = repmat (par_min,Nsmpl,1);
    else
        error (sprintf('Number of rows in par_min (%i) does not match number of samples (%i).',size(par_min,1),Nsmpl));
    end
end

% check par_max
if ~exist ('par_max','var')
    par_max = repmat (Inf,size(par0));
end
if isempty (par_max) % an empty vector was specified
    par_max = repmat (Inf,size(par0));
elseif size(par_max,2) ~= length(model_pars)
    error (sprintf('noblefit: number of columns in ''par_max'' (%i) must agree with the number of model parameters (%i) of the model function (%s).',size(par_max,2),length(model_pars),model_fun))
end
if size (par_max,1) ~= Nsmpl
    if size(par_max,1) == 1
        warning (sprintf('par_max contains only one row, but there are %i samples. Assuming same values of par_max for all samples...',Nsmpl));
        par_max = repmat (par_max,Nsmpl,1);
    else
        error (sprintf('Number of rows in par_max (%i) does not match number of samples (%i).',size(par_max,1),Nsmpl));
    end
end

% check if par_min < par_max
for i = 1:size(par_min,2)
    if any ( j = find (par_min > par_max) )
        error (sprintf('noblefit: min. allowed value for parameter %s is larger than max. allowed value in row %i.',model_pars{i},j(1)))
    end    
end

% check if par0 >= par_min and par0 <= par_max
for i = 1:size(par0,2)
    if any ( j = find (par0(:,i) < par_min(:,i)) ) % compare column by column
        error (sprintf('noblefit: par0 value is lower than par_min value for parameter %s in row %i.',model_pars{i},j(1)))
    end    
    if any ( j = find (par0(:,i) > par_max(:,i)) )
        error (sprintf('noblefit: par0 value is larger than par_max for parameter %s in row %i.',model_pars{i},j(1)))
    end    
end


% ****** PREPARE FITS ******

disp ('');
disp ('*********************************************************')
disp ('*** NOBLEFIT: PREPARING DATA AND MODEL FOR FITTING ***')
disp ('*********************************************************')


disp (sprintf('Fitting model function ''%s'' (file %s) using tracers %s.',model_fun,model_fun_file,trc));
model_args = ''; % string containing arguments for anonymous function call (used below)

% determine arguments for call to model function
i0 = iF = 1;
for i = 1:length(par_usage)
    switch par_usage(i)
        case 0 % use i-th parameter as fixed value
             model_args = sprintf ('%sP0(%i),',model_args,i0); % use the fixed value in par0
             i0 = i0 + 1;
        case 1
             model_args = sprintf ('%sPF(%i),',model_args,iF); % use the fit parameters (vector 'par')
             fitpar_names{iF} = model_pars{i};
             iF = iF + 1;
         otherwise
             error (sprintf('noblefit: unknown usage flag for model parameter %s: par_usage(%i) = %i.',model_pars{i},i,par_usage(i)));
     end % switch
end % for

% add tracers to model arguments
model_args = sprintf ('%s{',model_args); % open tracers cell string
for i = 1:Ntrcr
    model_args = sprintf ('%s''%s''',model_args,tracers{i});
    if i < Ntrcr % add comma
        model_args = sprintf ('%s,',model_args);
    else % close tracers cell string
        model_args = sprintf ('%s}',model_args);
    end
end

% determine how to translate the tracer_data struct into a matrix
tracer_data_map = repmat ([],1,Ntrcr);
for i = 1:Ntrcr
    k = find (strcmp(tracers{i},fieldnames(tracer_data)));
    if isempty(k)
        error (sprintf('noblefit: the tracer data do not contain the tracer ''%s''. Please check the names of the tracers included in the tracer data.',tracers{i}));
    end
    if length(k) > 1
        error (sprintf('noblefit: the tracer data contains more than one entry for tracer ''%s''.',tracers{i}));
    end
    tracer_data_map(i) = k; % the i-th tracer is in the k-th column of the data table
end

% check if data errors are non-zero (or even negative)
for i = 1:Nsmpl for j = 1:Ntrcr
    f = getfield (tracer_data,tracers{j});
    if f.err <= 0
        error (sprintf('noblefit: standard errors of observed data must be strictly positive numbers (sample %i, tracer ''%s''; : error = %g).',i,tracers{j}),f.err);
    end
end end

disp (''); fflush (stdout);


% ****** FIT THE MODEL ******

% ****** PREPARE FITS ******

disp ('');
disp ('**********************************')
disp ('*** NOBLEFIT: STARTING FITS ***')
disp ('**********************************')



% Initialize structs containing statistics info of ech fit:
chi2 = DF = pVal = repmat (NA,Nsmpl,1);

% move sample data from struct into matrices
X_val = X_err = repmat  (NA,Nsmpl,Ntrcr);
for i = 1:Nsmpl for j = 1:Ntrcr
    f = getfield (tracer_data(i),tracers{j}); % get struct value
    X_val(i,j) = f.val;
    X_err(i,j) = f.err;
end end

mdl = sprintf('%s(%s);',model_fun,model_args); % call to model function

par_val = par_err = repmat (NA,Nsmpl,NFpar);

for i = 1:Nsmpl % fit each sample
    disp ('');
    disp ('-----------------------------------'); fflush (stdout);
    disp ('');
    if isfield (tracer_data(i),'name')
        disp (sprintf('Sample %i of %i (%s):',i,Nsmpl,getfield(tracer_data(i)','name')));
    else
        disp (sprintf('Sample %i of %i:',i,Nsmpl));
    end
    fflush (stdout);

    % set up initial values etc. for i-th sample

    P0 = PF = PFnorm = PFmin = PFmax = [];
    for j = 1:length(par_usage)
        switch par_usage(j)
            case 0
                disp (sprintf('   %s = %g (fixed value)',model_pars{j},par0(i,j)));
                P0 = [P0,par0(i,j)]; % add fixed value to P0
            case 1
                disp (sprintf('   %s (fit parameter):',model_pars{j}))
                disp (sprintf('     - initial value\t= %g',par0(i,j)));
                disp (sprintf('     - allowed range\t= %g to %g',par_min(i,j),par_max(i,j)));

                PF    = [PF,par0(i,j)]; % add initial value to PF
                PFmin = [PFmin,par_min(i,j)]; % add min allowed value to PFmin
                PFmax = [PFmax,par_max(i,j)]; % add max allowed value to PFmin
    
                if isna(par_norm(i)) % guess a suitable value
                    PFnorm = [ PFnorm , par0(i,j) ];
                    if PFnorm(end) == 0
                        PFnorm(end) = 1;
                    end
                    disp (sprintf('     - Normalisation factor = %g (determined from initial value -- consider specifying values explicity!)',PFnorm(end)));
                else % use the specified value
                    PFnorm = [ PFnorm , par_norm(i,j) ];
                    disp (sprintf('     - Normalisation factor = %g',PFnorm(end)));
                end
        end % switch
    end % for
    fflush (stdout);

    k = find ( isfinite(X_val(i,:)) & isfinite(X_err(i,:)) );
    DF(i) = length(k) - NFpar;
    
    if DF(i) < 0 % there is no point in trying to fit anything
        disp (sprintf('   Not enough data to fit the model (DF = %i).',DF(i)));
        par_val(i,:) = par_err(i,:) = NA;
        chi2(i) = pVal(i) = NA;
        cov{i} = repmat (NA,NFpar,NFpar); % covariance matrix
        
    else
        
        
        % Define goal function that will be called by fminsearch to find the best-fit values of the model parameters:
        % - fminsearch starts with initial parameter values
        % - then varies parameter values to search for minimum. At the beginning, the size of variation will be similar for all fit parameters
        % Therefore: transform the values used by fminsearch before calling the model function with the actual model parameters (which may have very different values and ranges) like this:
        % (i) offset fit parameters to zero, i.e.: PF --> PF-PF0 (where PF is a vector of the fit parameters, and PF0 is the zero offset)
        % (ii) normalise the fit parameters so that their fminsearch variations are similar, i.e.: PF --> PF ./ PFnorm (where PFnorm is a vector of scaling factors)
        % (i)+(ii) transformed values used by fminsearch: PF = (PF-PF0)./PFnorm
        % ==> apply inverse transformation the parameter values used by fminsearch to parameter values for use with the model function: PF --> PF.*PFnorm + PF0
        % 
        % Define goal function for fminsearch as a call to the chi2 objective function with transformed parameter values:
        PF0 = PF; % initial values to start the fit, used as 'zero offset' in step (i) above
        goalfun = @(PF) nf_objfun (PF.*PFnorm+PF0,P0,PFmin,PFmax,mdl,X_val(i,:),X_err(i,:)); % call objective function with transformed fit parameters as described above        
        
        % turn off warnings about parameters that are outside the allowed range during fitting (this may happen too often if the fitter is trying to optimise outside the allowed range)
        warning ('off','nf_modelfun_params_out_of_bounds');
        warning ('off','atmos_gas_timerange');
        warning ('off','load-file-in-path'); % turn off warnings about loading data from atmospheric transient tracer data
        
        % do the fitting        
        [par_val(i,:),chi2(i)] = fminsearch(goalfun,repmat(0,size(PF0))); % find best-fit parameter values relative to PF0 (both offeset and normalised!)
        
        % turn warnings back on:
        warning ('on','nf_modelfun_params_out_of_bounds');
        warning ('on','atmos_gas_timerange');
        warning ('on','load-file-in-path'); % turn off warnings about loading data from atmospheric transient tracer data
        
        % evaluate fit results:
        par_val(i,:) = PF0 + par_val(i,:) .* PFnorm; % remove normalisation to obtain the right values
        pVal(i) = 1-chi2cdf (chi2(i),DF(i));
        
        % determine residuals ( MEAS - MODEL ) / (MEAS_ERR) at best-fit point (this will throw warnings that were turned off during fitting):
        PF = par_val(i,:);  eval (sprintf('M = %s;',mdl)); % evaluate model function at best-fit point
        res(i,:) = ( X_val(i,:) - M ) ./ X_err(i,:);
        
        % determine errors of the best-fit values using the Jacobian matrix of the model function (see Numerical Recipes, chapters 15.4 and 15.5)
        eval (sprintf('mfn = @(PF) %s',mdl)); % define model function with fit parameters (in PF) as variables
        A   = diag(1./X_err(i,:)) * jacobs(par_val(i,:),mfn); % design matrix of the fitting problem
        cov{i} = inv (A'*A); % this is the covariance matrix of the fit parameters. The diagonal contains the variances of the fit parameters
        par_err(i,:) = sqrt (diag(cov{i})); % standard uncertainties of the fit parameter values
        

    end % if DF < 0

    % show results:
    disp ('');
    disp ('   Statistics of fit:');
    disp (sprintf('   - chi2\t= %g',chi2(i)));
    disp (sprintf('   - DF\t\t= %g',DF(i)));
    disp (sprintf('   - p-value\t= %g',pVal(i)));

    disp ('');
    disp ('   Fit results:');
    for j = 1:NFpar
        disp (sprintf('   - %s\t\t= %g +/- %g',fitpar_names{j},par_val(i,j),par_err(i,j)));
    end
    disp ('');
    fflush (stdout);
        
end % for 'fit each sample'

disp (''); fflush (stdout);



% ****** SAY GOODBYE ******

disp ('*** NOBLEFIT: MODEL FITTING COMPLETED.')
disp ('');  fflush (stdout);




