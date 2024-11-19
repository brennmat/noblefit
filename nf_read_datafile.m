function data = nf_read_datafile (file,options);

% data = nf_read_datafile (file,options);
%
% Reads data from a formatted text file. The data needs to be organized in columns as follows:
% - Columns are assumed to be separated by tabs (or semicolons for *.csv files). Other delimiters may be specified using 'options'.
% - The first line must be a header line with names of the data in the columns.
% - The first column must contain sample names (treated as string).
% - The remaining columns must contain the data values (either numbers, NA, NaN, or empty).
% - If column titles include units (or anything else) in parentheses, the parentheses part is removed from the name (it may be useful to have the units in the data file, but the unit will be in the way for data formatting)
% - Data columns containing the data uncertainties (errors) are identified by adding 'err' somehere in the title. Example: if the Ne concentrations are given in column with title 'Ne', the column title of the corresponding errors could be 'Ne err', 'err. Ne', 'Ne_err', etc.
%
% INPUT:
% file: file name, may include path to file (string)
% options (optional): struct to provide options. May be useful to provide details about the format of the data file, may be useful to specify special file formats (e.g. ouptut from 4D database or input files for Franks noble fitter. Known options:
% - options.replace_zeros: if set, replace data values equal to zero by opt_replace_zeros (scalar)
% - options.filter_InvIsotopeRatios: if set to non-zero value (scalar), try to make sure that isotope ratios are such that the low-mass isotope is divided by the high-mass isotope (e.g., replace 40Ar/36Ar by 36Ar/40Ar).
%
% OUTPUT:
% data: array of structs with data values and corresponding errors. Every struct corresponds to one line in the data file. Fieldnames correspond to the column titles in the header line. Sample names are stored in field 'name'.
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

% set default values for options:
opt_replace_zeros           = [];
opt_filter_InvIsotopeRatios = 0; 
% check if options are specified, and parse them:
if exist ('options','var');
    f = fieldnames (options);
    for i = 1:length(f)
        switch upper(f{i})
            case 'REPLACE_ZEROS'
                opt_replace_zeros = getfield (options,f{i});
                
            case 'FILTER_INVISOTOPERATIOS'
                opt_filter_InvIsotopeRatios = getfield (options,f{i});
            
            otherwise
                warning (sprintf('Unknown options key: %s',f{i}));
        end % switch
    end %for
end

% set column delimiter:
[_,_,ext] = fileparts(file);
if toupper(ext) == '.CSV'
	delim = ';';
else
	delim = sprintf('\t');
end

% open text file for reading
[fid,msg] = fopen (file,'rt');
if fid < 0
    error (sprintf('Could not open text file for reading: %s',msg));
end

% read header and prepare data
H = '';
while length(H) == 0
	H = fgetl (fid);
	H = fliplr(deblank(fliplr(H))); % remove blanks at beginning of line
	if any ( k = findstr(H,'%') ) % remove comment (if any)
		k = k(1);
		if k > 1
			H = H(1:k);
		else
			H = '';
		end
	end
end
H = strsplit (H,delim);


Ncol = length (H); % number of columns

disp ('');
disp (sprintf('Parsing header (%i columns):',Ncol));
colinfo.field = '';
colinfo.is_err = 0;
colinfo = repmat (colinfo,Ncol,1);
ignorecols = [];
fields_invRatio = {};

for i = 1:Ncol % try to parse column titles

    n = H{i};
    disp (sprintf('Column %i:',i));
    disp (sprintf('   title: %s',n));
    
    if i == 1
        colinfo(i).field = 'name';

    else
        % prepare things
        N = upper (n);
        
        % parse parentheses
        u = [ findstr('(',n) findstr('[',n) findstr('{',n) findstr(')',n) findstr(']',n) findstr('}',n) ]; % find parentheses and stuff
        u = [ min(u) max(u) ]; % index to first and last parentheses
        if ~isempty(u) % remove parentheses part
            disp ('   removing parentheses...');
            n = [ n(1:u(1)-1) n(u(2)+1:end) ];
        end
        
        % parse hypens
        u = findstr ('-',n);
        n(u) = '_';
        
        % parse 'error' columns
        u = findstr('ERR',N);
        if ~isempty(u)
            colinfo(i).is_err = 1;
            u2 = [ findstr('ERR.',N) findstr('ERR_',N) findstr('_ERR',N) ];
            if isempty(u2)
                u = [ min(u) min(u)+3 ];
            else
                u = [ min(u2) min(u2)+4 ];
            end
            disp ('   This is an error column...');
            n = [ n(1:u(1)-1) n(u(2)+1:end)];
        end
        
        % remove spaces for further processing:
        u = ( n ~= ' ' );
        n = n(u);
        N = upper (n); % update N
        
        if findstr ('/',N) % parse known ratios
            disp ('   Parsing ratio to simple simple field name...')
            switch N
                case {'3HE/4HE','HE-3/HE-4'}
                    n = 'RHe';
                case {'4HE/3HE','HE-4/HE-3'}
                    if opt_filter_InvIsotopeRatios
                        n = 'RHe';
                        fields_invRatio{end+1} = 'RHe';
                    else
                        n = 'RHeInv';
                    end
                                       
                case {'20NE/22NE','NE-20/NE-22'}
                    n = 'RNe';
                case {'22NE/20NE','NE-22/NE-20'}
                    if opt_filter_InvIsotopeRatios
                        n = 'RNe';
                        fields_invRatio{end+1} = 'RNe';
                    else
                        n = 'RNeInv';
                    end
                case {'36AR/40AR','AR-36/AR-40'}
                    n = 'RAr';
                case {'40AR/36AR','AR-40/AR-36'}
                    if opt_filter_InvIsotopeRatios
                        n = 'RAr';
                        fields_invRatio{end+1} = 'RAr';
                    else
                        n = 'RArInv';
                    end
                case {'NE/XE'}
                    n = 'RNeXe';
                case {'AR/XE'}
                    n = 'RArXe';
                case {'KR/XE'}
                    n = 'RKrXe';

                otherwise
                    disp (sprintf('    Could not parse ratio (%s) to simple field name.',n));
            end % switch
        end % if
        
        % check if field name start with number, and parse them
        if ~( u = isnan(str2double(n(1)))) % field name starts with number, not good as field name in struct
            len = 1; % length of numeric string part
            while ~( u = isnan(str2double(n(len+1)))) % next string is also number
            	disp(len)
            	len += 1;
            end
            n = sprintf ('%s%s',n(1+len:end),n(1:len)); % move number from beginning to end
            N = upper (n); % update N
        end

        colinfo(i).field = deblank (n); % remove trailing blanks


    end % if i == 1

    
    % print column information:
    if i == 1
        disp (sprintf('   Corresponding field name in output: %s',colinfo(i).field))
    elseif colinfo(i).is_err
        disp (sprintf('   Corresponding field name in output: %s.err',colinfo(i).field))
    else
        disp (sprintf('   Corresponding field name in output: %s.val',colinfo(i).field))
    end
    fflush (stdout);
end % for 1:Ncol

fields_invRatio = unique(fields_invRatio);

disp ('');
disp ('Checking consistency of columns...')
for i = 2:Ncol
    if ~any( i == ignorecols )
        if ~colinfo(i).is_err % try to find matching column with error values
            found_errcol = 0;
            for j = 2:Ncol % search matching 'error' column
                if j ~= i
                    if strcmp (colinfo(j).field,colinfo(i).field) % there is another column with a matching field name
                        if ~colinfo(j).is_err % this is not an error column, so column j looks like a duplicate of column i
                            ignorecols = unique ([ ignorecols j ]);
                            disp (sprintf('   Column %i has the same field name as column %i. Ignoring data from column %i...',j,i,j));
                        else % this is a matching error column
                            found_errcol = 1;
                            break
                        end
                    end
                end            
            end % for
            if ~found_errcol
                ignorecols = unique ([ ignorecols i ]);
                disp (sprintf('   Column %i (field %s): could not find matching column with error values. Ignoring data from column %i...',i,colinfo(i).field,i));
            end
        end
    end % if ~any(find(i,ignorecols))
end



disp ('');

clear u
data = [];
L = fgetl (fid); LN = 2;
while L ~= -1 % while EOF is not reached
    disp (sprintf('Reading line %i...',LN));
    L = strsplit (L,delim);
    u.name = deblank(L{1}); % name of sample
    X = str2double({L{2:end}}); % data values of sample
    
    if ~isempty (opt_replace_zeros)
        X (find(X==0)) = opt_replace_zeros;
        disp (sprintf('Replaced zero values by %g...',opt_replace_zeros));
    end
    if any (k = isnan(X))
        X(k) = NA; % replace NaN by NA
        disp ('Replaced NaN values by NA...');
    end

    for i = 2:Ncol
        if ~any(i == ignorecols)
            if ~colinfo(i).is_err % value            
                eval (sprintf('u.%s.val = X(i-1);',colinfo(i).field))
            else % error
                eval (sprintf('u.%s.err = X(i-1);',colinfo(i).field))
            end            
        end % if ignorecol
    end % for
    
    data = [ data ; u ];

    fflush (stdout);
    
    L = fgetl (fid); LN += 1; % get next data line
end

if opt_filter_InvIsotopeRatios % deal with inverted isotope ratios
    if length(fields_invRatio) > 0
        for j = 1:length(fields_invRatio) % invert the values
            disp (sprintf('Inverting isotope ratio values (%s)...',fields_invRatio{j}))
            for i = 1:length(data)
                u = getfield (data(i),fields_invRatio{j});
                err = u.err / u.val;
                val = 1/u.val;
                err = err*val;
                eval (sprintf ('data(i).%s.val = val;',fields_invRatio{j}));
                eval (sprintf ('data(i).%s.err = err;',fields_invRatio{j}));
            end % for i = ...
        end % for j = ...
    end % if length(...)
end % if opt_filter...



disp ('Done.'); fflush (stdout);

fclose (fid); % close file after reading
