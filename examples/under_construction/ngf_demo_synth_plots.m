function ngf_demo_synth_plots (datafile)

% ngf_demo_synth_plots (datafile)
%
% Plots results from ngf_demo_synth_plots by loading a data file as produced by ngf_synth.m.
%
% EXAMPLE:
% ngf_demo_synth_plots ('ngf_synth_fitresults_232.mat')
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

load (datafile);

idx = find(iFit);
for i = 1:length(idx)
    figure; % prepare empty figure window
    __synth_plot (conc,parname{idx(i)},val(:,i),err(:,i),pVal);
end
end % function




% helper function:

function __synth_plot (u,f,y,y_err,pVal);
    x = [];
    for i = 1:length(u)
        eval (sprintf('x = [ x u(i).%s ];',f));
    end
%   X1 = min(x);
%   X2 = max(x);
%   line ( [X1 X2] , [X1 X2] );
%   hold on
    i = find (pVal >= 0.05); % acceptable fits 
    j = find (pVal < 0.05);  % poor fits
    
    if ~isempty(i)
        h1 = errorbar (x(i),y(i),y_err(i));
        set (h1,'linestyle','none','color','k','marker','o','markersize',3);
    end
    hold on
    
    if ~isempty(j)        
        h2 = errorbar (x(j),y(j),y_err(j));
        set (h2,'linestyle','none','color','r','marker','o','markersize',3); 
    end
    hold off
    
    axis ('equal');
    axis ('square');

    xlabel (sprintf('%s (original)',f));
    ylabel (sprintf('%s (fit)',f));
%    axis ( [ X1 - 0.05*(X2-X2) X1+0.05*(X2-X1) X1 - 0.05*(X2-X2) X1+0.05*(X2-X1) ] );
end
