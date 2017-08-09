%% qppOptions
%  Returns an option struct which can be passed to the presolver
%---------------------------------------------------------------------------------
%
% Returns a struct containing values for all options to be used within qpPresolver.
%
% Call
%    options = qppOptions( 'default' );
%    options = qppOptions( 'reliable' );
%    options = qppOptions( 'fast' );
% to obtain a set of default options or a pre-defined set of options tuned
% for reliable or fast QP preprocessing, respectively.
%
% Call
%    options = qppOptions( 'option1', value1, 'option2', value2, ... )
% to obtain a set of default options but with 'option1' set to value1 etc.
%
% Call
%    options = qppOptions( oldOptions, 'option1', value1, ... )
% to obtain a copy of the options struct oldOptions but with 'option1' set
% to value1 etc.
%
% Call
%    options = qppOptions( 'default',  'option1', value1, ... )
%    options = qppOptions( 'reliable', 'option1', value1, ... )
%    options = qppOptions( 'fast',     'option1', value1, ... )
% to obtain a set of default options or a pre-defined set of options tuned
% for reliable or fast QP preprocessing, respectively, but with 'option1' set to 
% value1 etc.
%
%
% qpPresolver features the following options:
% boundMode                     - Bounds that are passed to the solver:
%                                   'Medium': Medium bounds are passed to the solver
%                                   'Tightest': Tightest bounds are passed to the solver
% equalityTol                   - Tolerance > 0 when comparing two values for equality
% stabilityTol                  - Tolerance > 0 when performing sparsification (~ Gaussian
%                                   elimination) such that pivot element is not too small.
% feasibilityTol                - Tolerance > 0. Used for determination when the QP is to be
%                                   considered as infeasible.
% maxIter                       - Maximum number of iterations of the presolving loop.
% logfileLevel                  - Determines the amount of output written into the logfile.
%                                   0: No output.
%                                   1: Status of interface functions.
%                                   2: 1) + Status of presolving techniques.
%                                   3: 2) + Status of internal functions (not related to 
%                                      preprocessing techniques).
%                                   4: 3) + Status of internal functions (related to
%                                      preprocessing techniques).
%
% enableBoundTightening         - Enables (1) or disables (0) bound tightening of the
%                                      primal variables.
% enableDualConstraintsMethod   - Enables (1) or disables (0) dual constraints method.
% enableDuplicateColumnsMehtod  - Enables (1) or disables (0) duplicate columns method.
% enableEmptyColumnsMethod      - Enables (1) or disables (0) empty columns method.
% enablePrimalConstraintsMethod - Enables (1) or disables (0) the forcing, redundant and
%                                   infeasible primal constraints methods.
% enableScaling                 - Enables (1) or disables (0) scaling.
% enableSingletonColumnsMethod  - Enables (1) or disables (0) singleton columns method.
% enableSingletonRowsMethod     - Enables (1) or disables (0) singleton rows method.
% enableSparsificationMethod    - Enables (1) or disables (0) sparsification method.

function [options] = qppOptions(varargin)

firstIsStructOrScheme = 0;

if (nargin == 0)
    options = qppDefaultOptions();
else
    if (isstruct(varargin{1}))
        if (mod(nargin,2) ~= 1)
            error('qppOptions: Options must be specified in pairs!');
        end
        options = varargin{1};
        firstIsStructOrScheme = 1;
    else
        if (ischar(varargin{1}))
            if (mod(nargin,2) == 0)
                options = qppDefaultOptions();
            else
                if ( (nargin > 1) && (ischar(varargin(nargin))) )
                    error('qppOptions: Options must be specified in pairs!');
                end

                switch(lower(varargin{1}))
                    case 'default'
                        options = qppDefaultOptions();
                    case 'reliable'
                        options = qppReliableOptions();
                    case 'fast'
                        options = qppFastOptions();
                    otherwise
                        error(['qppOptions: Only the following option schemes are defined:', ...
                               '''default'', ''reliable'', ''fast''!']);
                end
                firstIsStructOrScheme = 1;
            end
        else
            error('qppOptions: First argument needs to be a string or an options struct!');
        end
    end
end

% Now set options based on user-defined values
for i = 1+firstIsStructOrScheme : 2 : nargin
    argName = varargin{i};
    argValue = varargin{i+1};

    if ( (isempty(argName)) || (~ischar(argName)) )
        error('qppOptions: Argument no. %d must be a non-empty string!', i);
    end
    
    if (~isfield(options, argName))
        error('qppOptions: Argument no. %d is an invalid option!', i);
    end
    
    if ( (strcmp(argName, 'boundMode') == false) && ...
            ( (ischar(argValue)) || (~isscalar(argValue)) ) )
        error('qppOptions: Argument no. %d must be a scalar!', i);
    end

    if (~ischar(argValue))
        eval(['options.', argName, ' = ', num2str(argValue), ';']);
    else
        eval(['options.', argName, ' = ''', argValue, ''';']);
    end
end

end



function [options] = qppDefaultOptions()

options = struct('boundMode',                      'Medium', ...
                 'equalityTol',                    1e-14, ...
                 'stabilityTol',                   1e-4, ...
                 'feasibilityTol',                 1e-8, ...
                 'maxIter',                        30, ...
                 'logfileLevel',                   2, ...
                 'enableBoundTightening',          1, ...
                 'enableDualConstraintsMethod',    1, ...
                 'enableDuplicateColumnsMethod',   1, ...
                 'enableEmptyColumnsMethod',       1, ...
                 'enablePrimalConstraintsMethod',  1, ...
                 'enableScaling',                  0, ...
                 'enableSingletonColumnsMethod',   1, ...
                 'enableSingletonRowsMethod',      1, ...
                 'enableSparsificationMethod',     1);
end


function [options] = qppReliableOptions()

options = qppDefaultOptions();

options.stabilityTol = 1e-2;
options.feasibilityTol = 1e-10;
options.enableSparsificationMethod = 0;

end


function [options] = qppFastOptions()

options = qppDefaultOptions();

options.enableDuplicateColumnsMethod = 0;
options.enableDualConstraintsMethod = 0;
options.enableBoundTightening = 0;
options.maxIter = 10;

end
