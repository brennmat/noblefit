function nf_save_fitresults (file,smpl,parnames,par_val,par_err,chi2,DF,p,fitcomment);

% function nf_save_fitresults (file,smpl,parnames,par_val,par_err,chi2,DF,p,fitcomment);
%
% Save fit results to ASCII file (e.g., for import to spreadsheet software).
%
% INPUT;
% T:        name of ASCII file (including file suffix)
% smpl:     vector of structs containing sample names (see also nf_read_datafile)
% parnames:	cell string containing names (and dimensions/units) of fit variables (used to determine header line)
% par_val:	matrix containing best-fit values of fit variables (as obtained from noblefit)
% par_err:	matrix containing errors of par_val (as obtained from noblefit)
% chi2, DF, p:		Vectors of Chi2, DF and p values (as obtained from noblefit)
% fitcomment: string with comments / notes about the fit(s) (optional)
%
% OUTPUT:
% (none)
%
% EXAMPLES:
% see files in Examples folder (e.g., EXAMPLE-2).
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

% check input
NS = length (smpl);
NP = length (parnames);

if size(par_val,2) ~= NP
	error (sprintf('Number of columns in par_val (%i) is different from length of parnames (%i).',size(par_val,2),NP));
end
if size(par_val,1) ~= NS
	error (sprintf('Number of rows in par_val (%i) is different from length of smpl (%i).',size(par_val,1),NS));
end
if size(par_err,2) ~= NP
	error (sprintf('Number of columns in par_err (%i) is different from length of parnames (%i).',size(par_err,2),NP));
end
if size(par_err,1) ~= NS
	error (sprintf('Number of rows in par_err (%i) is different from length of smpl (%i).',size(par_err,1),NS));
end
if length(chi2) ~= NS
	error (sprintf('Length of chi2 (%i) is different from length of smpl (%i).',length(chi2),NS));
end
if length(DF) ~= NS
	error (sprintf('Length of DF (%i) is different from length of smpl (%i).',length(DF),NS));
end
if length(p) ~= NS
	error (sprintf('Length of p (%i) is different from length of smpl (%i).',length(p),NS));
end


% open file for writing:
[fid,msg] = fopen (file,'wt');
if fid == -1
	error (sprintf('Could not open file for writing data (%s).',msg));
end

% write comment / info (if available)
if exist ('fitcomment','var')
	fprintf (fid,'COMMENT: %s\n\n',fitcomment);
end

% write header
fprintf (fid,'SAMPLE');
for i = 1:NP
	fprintf (fid,'\t%s\tErr. %s',parnames{i},parnames{i});
end
fprintf (fid,'\tChi2\tDF\tp');

% write data
for i = 1:NS
	fprintf (fid,'\n');
	fprintf (fid,smpl(i).name);
	for j = 1:NP
		fprintf (fid,'\t%g\t%g',par_val(i,j),par_err(i,j));
	end
	fprintf (fid,'\t%g\t%i\t%g',chi2(i),DF(i),p(i));
end

fclose (fid); % close file after reading