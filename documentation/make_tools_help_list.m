% this is an Octave m-file to make a LaTeX file containing the help messages of the NOBLEFIT m-files.

outfile = 'noblefit_manual_tools.tex';
out_fid = fopen(outfile,'wt');
pth = fileparts (which('noblefit'));
f = dir (sprintf('%s%s*.m',pth,filesep));

for i = 1:length(f)
    h = help (f(i).name);
    k = findstr (h,'DISCLAIMER');    
    if ~isempty(k)
        h = h(1:k-1);
    end
    k = findstr (h,'*********');    
    if ~isempty(k)
        h = h(1:k-1);
    end
    
    nam = f(i).name(1:end-2);
    nam2 = strrep(nam,'_','\_');
    
    fprintf (out_fid,'\\subsection{%s}\\label{sec:ref-%s}\n',nam2,nam);
    fprintf (out_fid,'\\begin{lstlisting}\n');
    fprintf (out_fid,'%s',h);
    fprintf (out_fid,'\n\\end{lstlisting}\n\n');

end % for i = ...
fclose (out_fid);