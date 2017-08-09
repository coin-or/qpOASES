%% install
%  Build qpPresolver library and install MEX interface.

function install(debugging)

isoctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);

% Check if 64 bit is supported
[~, maxInt] = computer;
is64bit = (maxInt > pow2(31));

% Check for debugging
if (nargin < 1)
    debugging = 0;
end


% ===================================================
% Compile library
% ===================================================

cmdline = 'make';

% Set compiler flags
flags = ['-Wall -Wextra -std=c99 -pedantic -Wmaybe-uninitialized ', ... 
         '-Wshadow -Wno-unused-function -finline-functions -fPIC '];
if (debugging)
    flags = [flags, '-g'];
else
    flags = [flags, '-O3'];
end

CFLAGS = [' CFLAGS="', flags, '"'];
cmdline = [cmdline, CFLAGS];

% Set Defines
defines = '';

% Set 64bit flag
if (is64bit)
    defines = [defines, '-DQPP_USE_INT64 '];
end

if (debugging)
    defines = [defines, '-DQPP_WRITE_LOGFILE '];
end

if (~isempty(defines))
    cmdline = [cmdline, ' DEFINES+="', defines, '"'];
end

disp(cmdline);

[stat, output] = system(cmdline);
if (stat ~= 0)
   fprintf('make failed:\n%s\n', output);
end


% ===================================================
% Compile MEX interface
% ===================================================

mexfiles = ['qpp.c ', ...
			'../../lib/libqpPresolver.a'];

% Mex command
mexcmd = 'mex';
if (is64bit)
	mexcmd = [mexcmd, ' -DQPP_USE_INT64'];
end

%	MATLAB does not allow C++ style comments (//) unless we use ISO C99 
%   (amongst other things); also pass additional flags to mex command.
if (~isoctave)
	mexcmd = [mexcmd, CFLAGS];
else
    mexcmd = [mexcmd, ' ', flags];
end

% Include directories
mexcmd = [mexcmd, ' -I../../include -I../../extern/mmio/include -I/usr/include/suitesparse'];

% Mex file and libraries
mexcmd = [mexcmd, ' ', mexfiles];

% Linker options
mexcmd = [mexcmd, ' -lrt -lsuitesparseconfig -lcholmod -lspqr -lumfpack'];
if (~isoctave)
    mexcmd = [mexcmd, ' -largeArrayDims -lmwblas -lmwlapack'];
end

eval(mexcmd);

end
