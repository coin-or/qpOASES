% example1
% Shows how to use the Matlab interface of qpPresolver to presolve LPs / QPs.

function [] = example1()

% QP can be found in example1.mat (Problem QBANDM from Maros-Meszaros test set)
load('example1.mat');


% Change options of the presolver (optionally). An ovierview of all options can be 
% found in qppOptions().
options = qppOptions('default', 'enableBoundTightening', 0, 'enableScaling', 1);

% Options can also be changed directly.
options.stabilityTol = 1e-4;


% Presolving QP; passing options is optional.
[id, qp, exitflag, iter, auxOutput] = qppPresolve(H, g, f, A, xl, xu, al, au, options);
if (exitflag ~= 0)
    error('QP could not be presolved. Error code %d', exitflag);
end

% Print information (time and dimension of QP)
fprintf('\nTime for presolving: %.6e seconds. Number of iterations: %d\n', ...
    auxOutput.time, iter);
fprintf('Dimension of presolved (original) QP:\n');
fprintf('\tNumber of variables        = %d (%d)\n', length(qp.g), length(g));
fprintf('\tNumber of lin. constraints = %d (%d)\n', size(qp.A,1), size(A,1));
fprintf('\tNumber of nonzeros in H    = %d (%d)\n', nnz(qp.H), nnz(H));
fprintf('\tNumber of nonzeros in A    = %d (%d)\n', nnz(qp.A), nnz(A));

% !!!!!!!!!!!!!!!
% Solve presolved QP (stored in struct qp) and retrieve primal-dual optimal 
% solution (px, py, pz)!
% !!!!!!!!!!!!!!!

% Here: "Dummy solution"
px = zeros(length(qp.g), 1);
py = zeros(size(qp.A, 1), 1);
pz = zeros(length(qp.g), 1);


% Compute primal-dual optimal solution (x, y, z) of ORIGINAL QP via postprocessing:
[exitflag, x, y, z] = qppPostsolve(id, px, py, pz);
if (exitflag ~= 0)
    error('Postprocessing falied. Error code %d', exitflag);
end

% If also optimal working sets of the presolved QP are given (for linear and bound
% constraints separately), then one can instead call
%
% [exitflag, x, y, z, wc, wb] = qppPostsolve(id, px, py, pz, pwc, pwb);
%
% where pwc is the optimal working set of the presolved QP corresponding to the linear
% constraints and pwb corresponding to the bound constraints.


% Do not forget to free memory allocated by the presolver entity:
qppFree(id);

end
