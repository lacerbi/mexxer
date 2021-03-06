function [lhs,rhs] = mexxer(filename,varargin)
%MEXXER Create template/interface of C code from MATLAB file.
%   MEXXER(FILENAME) generates C code from .m file FILENAME. MEXXER also 
%   creates a MATLAB script to test the MATLAB code vs MEX code.
%
%   MEXXER(FILENAME,'mex') only creates the MEX code.
%   MEXXER(FILENAME,'test') only creates the MATLAB test script.
%
%   MEXXER(...,1) overwrites the output file(s) if they already exist.

% Copyright (C) 2016 Luigi Acerbi
%
% This software is distributed under the GNU General Public License 
% (version 3 or later); please refer to the file LICENSE.txt, included with 
% the software, for details.
%
%   Author:     Luigi Acerbi
%   Email:      luigi.acerbi@gmail.com
%   Version:    24/Aug/2016

mexxerver = 'v0.2';

if nargin < 1 && nargout < 1
    fprintf('MEXXER version %s.\n\n', mexxerver);
    help mexxer.m; 
    return; 
end

% Last variable is overwrite
if nargin > 1 && isnumeric(varargin{nargin-1}) && isscalar(varargin{nargin-1})
    overwrite = varargin{nargin-1};
else
    overwrite = [];
end
if isempty(overwrite); overwrite = 0; end

writemexfun = [];
writetestfun = [];

if nargin > 1
    for i = 1:nargin-1
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'test'; writetestfun = 1; if isempty(writemexfun); writemexfun = 0; end
                case {'c','c++','mex'}; writemexfun = 1; if isempty(writetestfun); writetestfun = 0; end
                otherwise
                    error('mexxer:UnknownType',['Unknown file type ''' varargin{i} '''.']);
            end
        end
    end
end

if isempty(writemexfun); writemexfun = 1; end
if isempty(writetestfun); writetestfun = 1; end

% Ensure that all memory is clean
clear functions; fclose('all');

% Open input file
[filepath,name,ext] = fileparts(filename);
if isempty(ext); ext = '.m'; end
filename = fullfile(filepath,[name,ext]);
if ~exist(filename,'file')
    error(['File ' filename ' does not exist.']);
end

% Read function from input file
fin = fopen(filename,'r');
while 1        
    fundef = fgetl(fin);
    if ~isemptyline(fundef); break; end
end

% Parse function definition
fundef = regexprep(fundef,'[,()\[\]]',' ');
fundef = regexprep(fundef,' +',' ');    % Trim multiple whitespace

idx = find(fundef == '=');
if ~isempty(idx)
    lhs_list = strread(fundef(1:idx-1),'%s')';
    lhs_list(1) = []; % Remove 'function' token
    rhs_list = strread(fundef(idx+1:end),'%s')';
else
    lhs_list = [];
    rhs_list = strread(fundef,'%s')';
end    

% Store function description
while 1
    desc{1} = fgetl(fin);
    if ~isemptyline(desc{1}); break; end
end
idx = find(desc{end} == '%',1);
desc{end}(1:idx) = [];
while 1
    desc{end+1} = fgetl(fin);
    if isemptyline(desc{end}) || all(desc{end} ~= '%'); break; end
    idx = find(desc{end} == '%',1);
    desc{end}(1:idx) = [];
end
desc(end) = []; % Discard last line
fclose(fin);    % Close file

nlhs = numel(lhs_list);
funcname = rhs_list{1};
rhs_list(1) = [];
nrhs = numel(rhs_list);

% Parse arguments information from file
[rhssize,rhstype,rhsdesc] = parsevariables(desc,'input',nrhs);
[lhssize,lhstype,lhsdesc] = parsevariables(desc,'output',nlhs);

% Fill in argument struct
rhs = buildargstruct(rhs_list,rhssize,rhstype,rhsdesc,'input');
lhs = buildargstruct(lhs_list,lhssize,lhstype,lhsdesc,'output');
allargs = [lhs,rhs];

% Find free variables (input array sizes)
inputvars = {rhs.name};
freevars = [];
for i = 1:nrhs
    n = numel(rhs(i).sizes);
    rhs(i).getsize = zeros(1,n);
    if n == 1; continue; end     % Scalar, nothing interesting
    for j = 1:n
        s = rhs(i).sizes{j};
        if isfinite(str2double(s)); continue; end   % Constant
        rhs(i).vars{j} = strread(regexprep(s,'[,()\[\] *+/=-]',' '),'%s');
        for k = 1:numel(rhs(i).vars{j})
            if any(strcmp(inputvars,rhs(i).vars{j}{k})); continue; end
            freevars{end+1} = rhs(i).vars{j}{k};
            if numel(rhs(i).vars{j}) == 1; rhs(i).getsize(k) = 1; end 
        end
    end
end
freevars = unique(freevars);
freevars_defined = zeros(1,numel(freevars));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START WRITING ON FILE

if writemexfun

    if ~strcmp(name,funcname)
        warning(['Function file name ''' name ''' and function name ''' funcname ''' do not match.']);
    end

    fileout = fullfile(filepath,[name,'_mex.c']);
    if exist(fileout,'file') && ~overwrite
        error([fileout ' already exists. Call MEXXER(FILENAME,1) if you want to overwrite an existing file (be careful).']);
    end
    fout = fopen(fileout,'w');

    % File header
    fprintf(fout, '#include "mex.h"\n#include "math.h"\n#include "matrix.h"\n#include "float.h"\n\n');

    % Function definition
    fprintf(fout, ['/*\n * ' [name,'_mex.c'] '\n *\n']);
    for i = 1:numel(desc)
        fprintf(fout, ' *%s\n', desc{i});
    end
    % fprintf(fout, [' *\n */\n\n']);

    % Timestamp and info
    fprintf(fout, ' *\n * This is a MEX-file for MATLAB.\n * Template C code generated on %s with MEXXER %s \n * (https://github.com/lacerbi/mexxer).\n */\n\n', date, mexxerver);

    % Macros and definitions
    fprintf(fout, '/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */\n#define ARGSCHECK 1\n\n');

    % Specific function
    fprintf(fout, 'void %s( ',funcname);
    for i = 1:nlhs
        if i > 1; prefix = ', '; else prefix = ''; end
        fprintf(fout, '%s%s %s%s', prefix, lhs(i).type, repmat('*',[1,lhs(i).pointer]), lhs(i).name);
    end
    for i = 1:nrhs
        if nlhs > 0 || i > 1; prefix = ', '; else prefix = ''; end
        fprintf(fout, '%s%s %s%s', prefix, rhs(i).type, repmat('*',[1,rhs(i).pointer]), rhs(i).name);
    end
    for i = 1:numel(freevars)
        if nlhs > 0 || nrhs > 0 || i > 1; prefix = ', '; else prefix = ''; end
        fprintf(fout, '%smwSize %s', prefix, freevars{i});    
    end
    fprintf(fout,' )\n{\n\t\n\t/* Write your main calculations here... */\n\t\n}\n\n');

    % Main function
    fprintf(fout, '/* the gateway function */\n');

    fprintf(fout, 'void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )\n{\n');

    vartypes = unique({allargs.type});
    for i = 1:numel(vartypes)    
        fprintf(fout,'\t%s',vartypes{i});
        first = 1;
        for j = 1:numel(allargs)
            if strcmp(allargs(j).type,vartypes{i});
                pointer = repmat('*',[1,allargs(j).pointer]);
                if first; first = 0; prefix = ' '; else prefix = ', '; end
                fprintf(fout,'%s%s%s',prefix,pointer,allargs(j).name);
            end
        end
        fprintf(fout,';\n');    
    end

    % Add vector arrays for input and output arrays sizes
    fprintf(fout,'#if ( ARGSCHECK==0 )\n');
    first = 1;
    for i = 1:nrhs
        n = numel(rhs(i).sizes);
        if n > 1 && any(rhs(i).getsize)
            if first; fprintf(fout, '\tmwSize'); first = 0; prefix = ' '; else prefix = ', '; end
            fprintf(fout, '%s*dims_%s', prefix, rhs(i).name);
        end
    end
    if first == 0; fprintf(fout,';\n'); end
    fprintf(fout,'#else /* ( ARGSCHECK!=0 ) */ \n');
    first = 1;
    for i = 1:nrhs
        n = numel(rhs(i).sizes);
        if n > 1
            if first; fprintf(fout, '\tmwSize'); first = 0; prefix = ' '; else prefix = ', '; end
            fprintf(fout, '%s*dims_%s', prefix, rhs(i).name);
        end
    end
    if first == 0; fprintf(fout,';\n'); end
    fprintf(fout,'#endif /* ( ARGSCHECK!=0 ) */ \n');

    first = 1;
    for i = 1:nlhs
        n = numel(lhs(i).sizes);
        if n > 2
            if first; fprintf(fout, '\tmwSize'); first = 0; prefix = ' '; else prefix = ', '; end
            fprintf(fout, '%sdims_%s[%d]', prefix, lhs(i).name, n);
        end
    end
    if first == 0; fprintf(fout,';\n'); end

    % Free variables (matrix sizes)
    for i = 1:numel(freevars)
        if i == 1; fprintf(fout, '\tmwSize'); prefix = ' '; else prefix = ', '; end
        fprintf(fout,'%s%s',prefix,freevars{i});
        if i == numel(freevars); fprintf(fout, ';\n'); end
    end
    fprintf(fout,'\n');

    % Check for number of arguments (this is always done)
    fprintf(fout,'\t/*  check for proper number of arguments */\n\t/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt\n\t   within an if statement, because it will never get to the else\n\t   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks\n\t   you out of the MEX-file) */\n');
    fprintf(fout,'\tif ( nrhs<%d || nrhs>%d )\n', nrhs, nrhs); 
    fprintf(fout,'\t\tmexErrMsgIdAndTxt( "MATLAB:%s:invalidNumInputs",\n\t\t\t"%s inputs required.");\n', funcname, cardinal(nrhs,1));
    fprintf(fout,'\tif ( nlhs<%d || nlhs>%d )\n', nlhs, nlhs); 
    fprintf(fout,'\t\tmexErrMsgIdAndTxt( "MATLAB:%s:invalidNumOutputs",\n\t\t\t"%s outputs required.");\n\n', funcname, cardinal(nlhs,1));

    % Get inputs
    for i = 1:nrhs
        argdescription(rhs(i),i,'input',fout);
        if strcmp(rhs(i).sizes{1},'scalar')
            fprintf(fout,'\t%s = (%s) mxGetScalar(prhs[%d]);\n\n', rhs(i).name, rhs(i).fulltype, i-1);
        else
            fprintf(fout,'\t%s = (%s) mxGetPr(prhs[%d]);\n', rhs(i).name, rhs(i).fulltype, i-1);
            if any(rhs(i).getsize)
                fprintf(fout, '\tdims_%s = (mwSize*) mxGetDimensions(prhs[%d]);\n', rhs(i).name, i-1);
                % Check if any free variable needs to be read here
                for j = 1:numel(rhs(i).vars)
                    if numel(rhs(i).vars{j}) > 1; continue; end
                    idx = find(strcmp(rhs(i).vars{j},freevars),1);
                    if isempty(idx) || freevars_defined(idx); continue; end
                    fprintf(fout, '\t%s = dims_%s[%d];\n', freevars{idx}, rhs(i).name, j-1);
                    freevars_defined(idx) = 1;
                end
            end
            fprintf(fout,'\n');
        end
    end

    % Do input argument checking
    fprintf(fout,'\t/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */\n#if ( ARGSCHECK==1 )\n');
    first = 1;
    for i = 1:nrhs
        if strcmp(rhs(i).sizes{1},'scalar')
            if first; first = 0; else fprintf(fout,'\n'); end
            fprintf(fout,'\t\tif ( !mxIsDouble(prhs[%d]) || mxIsComplex(prhs[%d]) || (mxGetN(prhs[%d])*mxGetM(prhs[%d])!=1) )\n', i-1, i-1, i-1, i-1);
            fprintf(fout,'\t\t\tmexErrMsgIdAndTxt("MATLAB:%s:%sNotScalar", "Input %s must be a scalar.");\n', funcname, rhs(i).name, upper(rhs(i).name));
        else
            if first; first = 0; else fprintf(fout,'\n'); end
            if all(rhs(i).getsize == 0)
                fprintf(fout, '\t\tdims_%s = (mwSize*) mxGetDimensions(prhs[%d]);\n', rhs(i).name, i-1);
            end
            fprintf(fout,'\t\tif ( !mxIsDouble(prhs[%d]) || mxIsComplex(prhs[%d]) )\n\t\t\t\tmexErrMsgIdAndTxt("MATLAB:%s:%sNotReal", "Input %s must be real.");\n', i-1, i-1, funcname, rhs(i).name, upper(rhs(i).name));
            for k = 1:numel(rhs(i).sizes)
                fprintf(fout,'\t\tif ( dims_%s[%d] != ((mwSize) (%s)) )\n\t\t\tmexErrMsgIdAndTxt("MATLAB:%s:%sWrongSize", "The %s dimension of input %s has the wrong size (should be %s).");\n', rhs(i).name, k-1, rhs(i).sizes{k}, funcname, rhs(i).name, ordinal(k), upper(rhs(i).name), rhs(i).sizes{k});
            end        
        end
    end
    fprintf(fout,'#endif /* ( ARGSCHECK==1 ) */ \n\n');

    % Prepare outputs and pointers
    for i = 1:nlhs
        argdescription(lhs(i),i,'output',fout);
        n = numel(lhs(i).sizes);
        switch n
            case 1
                fprintf(fout,'\tplhs[%d] = mxCreateDoubleScalar(0.);\n', i-1);            
            case 2
                fprintf(fout,'\tplhs[%d] = mxCreateDoubleMatrix((mwSize) (%s), (mwSize) (%s), mxREAL);\n', i-1, lhs(i).sizes{1}, lhs(i).sizes{2});
            otherwise
                for j = 1:n; fprintf(fout, '\tdims_%s[%d] = (mwSize) (%s);\n', lhs(i).name, j-1, lhs(i).sizes{j}); end
                fprintf(fout,'\tplhs[%d] = mxCreateNumericArray(%d, dims_%s, mxDOUBLE_CLASS, mxREAL);\n', i-1, n, lhs(i).name);
        end
        fprintf(fout,'\t%s = mxGetPr(plhs[%d]);\n\n', lhs(i).name, i-1);
    end

    % Call subroutine
    fprintf(fout, '\t/* Call the C subroutine */\n\t%s(', funcname);
    for i = 1:nlhs
        if i > 1; prefix = ', '; else prefix = ''; end
        fprintf(fout, '%s%s', prefix, lhs(i).name);
    end
    for i = 1:nrhs
        if nlhs > 0 || i > 1; prefix = ', '; else prefix = ''; end
        fprintf(fout, '%s%s', prefix, rhs(i).name);
    end
    for i = 1:numel(freevars)
        if nlhs > 0 || nrhs > 0 || i > 1; prefix = ', '; else prefix = ''; end
        fprintf(fout, '%s%s', prefix, freevars{i});    
    end
    fprintf(fout, ');\n\n');

    fprintf(fout,'}\n');

    fclose(fout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST FUNCTION

if writetestfun

    fileout = fullfile(filepath,[name,'_test.m']);
    if exist(fileout,'file') && ~overwrite
        error([fileout ' already exists. Call MEXXER(FILENAME,1) if you want to overwrite an existing file (be careful).']);
    end
    fout = fopen(fileout,'w');

    % Timestamp and info
    fprintf(fout, '%%%s\n%%   Test script for MEX-file %s.\n%%\n', upper([name '_test']), name);
    fprintf(fout, '%%   Template MATLAB code generated on %s with MEXXER %s \n%%   (https://github.com/lacerbi/mexxer).\n\n', date, mexxerver);
    
    % Default parameters
    fprintf(fout, 'TolErr = sqrt(eps);\t%% Maximum error tolerance per array element\n\n');
    
    % Initialize sizes
    fprintf(fout, '%% Define array sizes (choose reasonable values)\n');
    vars = [];
    for i = 1:numel(rhs)
        if ~isempty(rhs(i).vars); vars = [vars, rhs(i).vars{:}]; end
    end
    vars = unique(vars);
    for i = 1:numel(vars)
        fprintf(fout, '%s = %d;\n', vars{i}, 50 + 50*i);
    end

    % Randomly initialize variables
    fprintf(fout, '\n%% Randomly initialize input variables\n%% (or write here alternative initializations)\n');
    for i = 1:numel(rhs)
        % Skip variables if already initialized
        if any(strcmp(rhs(i).name,vars)); continue; end

        switch rhs(i).type
            % case 'double'
            case 'int'; rngmethod = 'randi(10,';
            otherwise; rngmethod = '10*rand(';
        end
        fprintf(fout, '%s = %s[', rhs(i).name, rngmethod);
        if numel(rhs(i).sizes) == 1
            fprintf(fout,'1,1');
        else    
            fprintf(fout,'%s',rhs(i).sizes{1});
            for j = 2:numel(rhs(i).sizes)
                fprintf(fout,',%s',rhs(i).sizes{j});
            end
        end
        fprintf(fout,']);\t%%%s\n',rhs(i).desc);  
    end

    basestring = ['Testing ' name ':'];
    fprintf(fout,'\nfprintf(''%s\\n'');\n',repmat('=',[1,numel(basestring)]));
    fprintf(fout,'fprintf(''%s\\n'');\n', basestring);
    fprintf(fout,'fprintf(''%s\\n'');\n',repmat('=',[1,numel(basestring)]));

    fprintf(fout,'\n%% Call MATLAB and MEX functions\n');

    % Call MATLAB function
    fprintf(fout,'tic; ');
    if numel(lhs) > 0
        fprintf(fout,'[%s', lhs(1).name);
        for i = 2:numel(lhs); fprintf(fout,',%s', lhs(i).name); end
        fprintf(fout,'] = ');
    end
    fprintf(fout,'%s(',name);
    if numel(rhs) > 0
        fprintf(fout,'%s', rhs(1).name);
        for i = 2:numel(rhs); fprintf(fout,',%s', rhs(i).name); end
    end
    fprintf(fout,'); t = toc;\n');

    % Call MEX function
    fprintf(fout,'tic; ');
    if numel(lhs) > 0
        fprintf(fout,'[%s_mex', lhs(1).name);
        for i = 2:numel(lhs); fprintf(fout,',%s_mex', lhs(i).name); end
        fprintf(fout,'] = ');
    end
    fprintf(fout,'%s_mex(',name);
    if numel(rhs) > 0
        fprintf(fout,'%s', rhs(1).name);
        for i = 2:numel(rhs); fprintf(fout,',%s', rhs(i).name); end
    end
    fprintf(fout,'); t_mex = toc;\n\n');

    fprintf(fout,'%% Correctness check\n');
    for i = 1:numel(lhs)
        fprintf(fout, '%s_err = sum(abs(%s(:)-%s_mex(:)));\n', lhs(i).name, lhs(i).name, lhs(i).name);
        fprintf(fout, ...
            'fprintf(''Total error (%s): %%g\\n'', %s_err);\n', ...
            lhs(i).name, lhs(i).name);
        fprintf(fout,'if %s_err > TolErr*numel(%s);\n\terror(''mexxer:tooLargeError'',''Correctness check failed. Error too large in %s.'');\nend\n', lhs(i).name, lhs(i).name, lhs(i).name);
    end

    fprintf(fout,'\n%% Runtime analysis\n');
    fprintf(fout,'fprintf(''Time for MATLAB code: %%.3f s\\n'', t);\n');
    fprintf(fout,'fprintf(''Time for MEX file: %%.3f s\\n'', t_mex);\n');
    fprintf(fout,'fprintf(''Speed gain: %%.2f\\n'', t/t_mex);\n');


    fclose(fout);
end

end

%--------------------------------------------------------------------------
function tf = isemptyline(line)
%ISEMPTYLINE Return TRUE if line is empty space.
tf = (isnumeric(line) && line == -1) | isempty(regexprep(line,'[\b\f\n\r\t ]',''));
end

%--------------------------------------------------------------------------
function [tf,idx] = issepline(line)
%ISSEPLINE Return TRUE if line is a separator (also return start of sep).
idx = min([Inf,strfind(line,'==='),strfind(line,'%%%')]);
tf = isfinite(idx);
end

%--------------------------------------------------------------------------
function [varsize,vartype,vardesc] = parsevariables(desc,section,nvars)
%PARSEVARIABLES Parse argument description for size and type.

varsize = [];
vartype = [];
vardesc = [];

% Parse variables
found = 0;
for i = 1:numel(desc)
    line = desc{i};
    if ~issepline(line); continue; end
    idx = min([Inf,strfind(lower(line),section)]);
    if isfinite(idx); found = 1; break; end    
end

ivars = 0;
if found
    for i = i+1:numel(desc)
        line = desc{i};
        if issepline(line); break; end
        line = regexprep(line,'[%]',' ');
        if isemptyline(line); break; end
        
        ivars = ivars + 1;
        if ivars > nvars
            error(['Too many ' section ' arguments in file description. Each argument should have one line.']);
        end
        
        % Find sizes, last square brackets
        idx1 = find(line == '[',1,'last')+1;
        idx2 = idx1 - 1 + find(line(idx1+1:end) == ']',1);     
        if isempty(idx1) || isempty(idx2)
            error(['Cannot find size information for ' section ' #' num2str(ivars) ' in ' filein '.']);
        end        
        varsize{ivars} = line(idx1:idx2);
        
        % Find type, last round brackets
        idx1 = find(line == '(',1,'last')+1;
        idx2 = idx1 - 1 + find(line(idx1+1:end) == ')',1);        
        if isempty(idx1) || isempty(idx2)
            error(['Cannot find type information for ' section ' #' num2str(ivars) ' in ' filein '.']);
        end        
        vartype{ivars} = line(idx1:idx2);
        
        % Find description, until last round/square bracket
        idx2 = min(find(line == '(',1,'last'),find(line == '[',1,'last'))-1;
        vardesc{ivars} = line(1:idx2);
    end
end

if ivars < nvars
    error('Too few input arguments in file description. Each argument should have one line.');
end

end

%--------------------------------------------------------------------------
function astruct = buildargstruct(alist,asize,atype,adesc,section)
%BUILDARGSTRUCT Build argument structures after parsing from file.

n = numel(alist);
for i = 1:n
    % Remove whitespace from name (if any)
    astruct(i).name = regexprep(alist{i},'[\b\f\n\r\t ]',' ');
    
    % Scan argument sizes
    temp = regexprep(asize{i},'[\b\f\n\r\t\[\],; ]',' ');
    while 1
        idx = unique([strfind(temp,'-by-'),strfind(temp,'-BY-')]);
        if isempty(idx); break; end
        temp(idx(1)) = ' ';
        temp(idx(1)+1:idx(1)+3) = [];
    end
        
    errstr = ['Cannot parse ' section ' argument #' num2str(i) ' size in function description.'];
    if isempty(temp); error(errstr); end
    sizes = strread(temp,'%s')';
    if isempty(sizes); error(errstr); end
    
    if numel(sizes) == 1
        if any(strcmpi(sizes{1},{'1','scalar'}))
            astruct(i).sizes{1} = 'scalar';
        else
            error([upper(section(1)) section(2:end) ' argument #' num2str(i) ' size has only one dimension. Write ''1'' or ''scalar'' for a scalar value.']);
        end
    else
        if all(strcmp(sizes,'1'))
            astruct(i).sizes{1} = 'scalar';
        else            
            astruct(i).sizes = sizes;
        end
    end
    
    % Scan argument type
    temp = atype{i};
    errstr = ['Cannot parse ' section ' argument #' num2str(i) ' type in function description.'];
    if isempty(temp); error(errstr); end
    
    % Check if it is a pointer
    astruct(i).pointer = max(sum(temp == '*'), numel(astruct(i).sizes) > 1);

    if astruct(i).pointer && strcmp(astruct(i).sizes{1},'scalar')
        error([upper(section(1)) section(2:end) ' argument #' num2str(i) ' is a scalar pointer. Just declare it as a scalar.']);
    end
    
    temp = regexprep(temp,'[\b\f\n\r\t\[\],;* ]','');
    
    switch lower(temp)
        case {'int','integer'}; astruct(i).type = 'int';
        case 'double'; astruct(i).type = 'double';
        case 'float'; astruct(i).type = 'float';
        case {'size','mwsize'}; astruct(i).type = 'mwSize';
        case {'index','signedindex','mwsignedindex'}; astruct(i).type = 'mwSignedIndex';
        otherwise
            warnstr = ['Unknown variable type ''' temp ''' for ' section ' argument #' num2str(i) ' in function description.'];
            warning(warnstr);
            astruct(i).type = temp;
    end
    
    if astruct(i).pointer && ~strcmp(astruct(i).type,'double')
        warning([upper(section(1)) section(2:end) ' argument #' num2str(i) ' is a non-double pointer. This might cause problems.']);
    end
    astruct(i).fulltype = [astruct(i).type repmat('*',[1,astruct(i).pointer])];
    
    astruct(i).desc = strtrim(adesc{i});
end

end
%--------------------------------------------------------------------------
function argdescription(arg,n,section,fout)
%ARGDESCRIPTION Write argument description

switch lower(section)
    case 'input'
        fprintf(fout,'\t/* Get %s input (', ordinal(n));
    case 'output'
        fprintf(fout,'\t/* Pointer to %s output (', ordinal(n));
end
fprintf(fout,'%s, ', upper(arg.name));

for j = 1:numel(arg.sizes)
    if j > 1; prefix = '-by-'; else prefix = ''; end
    fprintf(fout,'%s%s',prefix,arg.sizes{j});
end
fprintf(fout,' %s) */\n', arg.type);

end

%--------------------------------------------------------------------------
function s = cardinal(n,firstupper)
%CARDINAL Return a string with a cardinal number
if nargin < 2 || isempty(firstupper); firstupper = 0; end

switch n
    case 0; s = 'zero';
    case 1; s = 'one';
    case 2; s = 'two';
    case 3; s = 'three';
    case 4; s = 'four';
    case 5; s = 'five';
    case 6; s = 'six';
    case 7; s = 'seven';
    case 8; s = 'eight';
    case 9; s = 'nine';
    case 10; s = 'ten';
    otherwise; s = num2str(n);
end
if firstupper; s(1) = upper(s(1)); end

end

%--------------------------------------------------------------------------
function s = ordinal(n,firstupper)
%ORDINAL Return a string with an ordinal number
if nargin < 2 || isempty(firstupper); firstupper = 0; end

%ordinal = {'1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th', ...
%    '11th','12th','13th','14th','15th','16th','17th','18th','19th','20th'};

switch n
    case 0; s = 'zeroth';
    case 1; s = 'first';
    case 2; s = 'second';
    case 3; s = 'third';
    case 4; s = 'fourth';
    case 5; s = 'fifth';
    case 6; s = 'sixth';
    case 7; s = 'seventh';
    case 8; s = 'eighth';
    case 9; s = 'ninth';
    case 10; s = 'tenth';
    otherwise; s = [num2str(n) 'th'];
end
if firstupper; s(1) = upper(s(1)); end

end