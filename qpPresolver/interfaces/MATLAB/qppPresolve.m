%% qppPresolve
% Presolves the given QP
%---------------------------------------------------------------------------------
%
% A QP of the form
%
%       min  1/2 x' * H * x + g' * x + f
%       s.t. al <= A*x <= au
%            xl <= x <= xu
%
% will be preprocessed (several transformations are applied to decrease the size of the
% QP (amongst other things).
%
% Call
%   [id, qp, exitflag, iter, auxOutput] = qppPresolve(H, g, f, A, xl, xu, al, au{, options})
%
% to presolve the above mentioned QP. H must be a symmetric matrix and all vectors g, xl,
% xu, al and au must be passed as column vectors. The gradient vector g must not be empty,
% but all other quantities may be empty, e.g. A = [], al = [] and au = [] if no linear
% constraints are present. 
% If empty bounds are passed to the presolver, e.g. xl = [] or al = [], then these 
% bounds will considered as -infinity. Analogously, if xu = [] or au = [], then these
% bounds are considered as +infinity.
% Options can be generated using the qppOptions() function, otherwise default values 
% are used. Matrices A and H will be converted into sparse format if they are not
% already sparse.
%
% Outputs (only id is mandatory):
%   id              - Unique ID for the newly created presolver entity. Must be known for
%                     later postprocessing, etc.
%   qp              - Struct containing the presolved QP (i.e. H, g, f, ...).
%                     Access via qp.H, qp.g, ...
%   exitflag        - Error code (for further details see errorcodes.h in
%                     include/qpPresolver). On success, exitflag = 0.
%   iter            - Number of iterations of the presolving loop.
%   auxOutput       - Struct containing auxiliary outputs as described below.
%
% The auxOutput struct contains the following entries:
%   time            - Wall clock time needed for preprocessing the QP.

function [id, qp, exitflag, iter, auxOutput] = qppPresolve(H, g, f, A, xl, xu, al, au, options)

if (nargin < 8)
    error('qppPresolve: Not enough input parameters!');
end

if (nargin < 9)
    options = [];
end
qp = [];
auxOutput = struct('time', 0.0);

if (~issparse(A))
    A = sparse(A);
end

if (~issparse(H))
    H = sparse(H);
end

% We only need the lower triangular part of the Hessian, but if the user
% passes a "full" Hessian matrix (with upper triangular part), then also a
% "full" matrix will be returned.
nnzH = nnz(H);
H = tril(H);
returnFullHessian = (nnz(H) < nnzH);

% Get dimension of the QP
[m,~] = size(A);
n = length(g);

% Create new qpPresolver entity
[id, exitflag] = qpp(0, m, n, nnz(A), nnz(H), 1000);
if (exitflag ~= 0)
    warning('Could not create new qpPresolver entity!');
    return;
end

% Initialize the newly created entity with given QP
exitflag = qpp(1, id, A, H, g, f, xl, xu, al, au);
if (exitflag ~= 0)
    warning('Could not initialize qpPresolver entity with given QP!');
    return;
end

% Set user options (if given) before presolving
if (~isempty(options))
    qpp(30, id, options);
end

% Now presolve the QP, measure wall clock time
tic();
[exitflag, iter] = qpp(10, id);
presolveTime = toc();
if (exitflag ~= 0)
    warning('Presolving unsuccessful, cf. exitflag!');
    return;
end

% Retrieve the presolved QP
[exitflag, m, n, Airn, Ajcn, Ax, Hirn, Hjcn, Hx, g, f, xl, xu, al, au] = ...
    qpp(15, id);
if (exitflag ~= 0)
    warning('Could not retrieve presolved QP!');
    return;
end

% Converting 0-based coordinate scheme to 1-based CSC scheme!
A = sparse(double(Airn)+1, double(Ajcn)+1, Ax, m, n);
H = sparse(double(Hirn)+1, double(Hjcn)+1, Hx, n, n);

if (returnFullHessian)
    H = H + transpose(tril(H,-1));
end

qp = struct('H', H, 'g', g, 'f', f, 'A', A, 'xl', xl, 'xu', xu, 'al', al, 'au', au);
auxOutput.time = presolveTime;

end