%% qppPostsolve
%  Perform postsolving operations based on the presolved QP
%---------------------------------------------------------------------------------
%
% Returns a primal-dual optimal solution (x, y, z) of the original QP based on a
% primal-dual optimal solution of the presolved QP (px, py, pz). Note that y and py refer
% to the multipliers corresponding to the linear constraints and z and pz corresponding to
% the bound constraints.
% Furthermore, the postprocessing routine optionally tries to determine an optimal working
% set (wc, wb) of the original QP based on an optimal working set of the presolved QP 
% (pwc, pwb). Here, wc and pwc refer to the working set corresponding to the linear
% constraints and wb and pwb corresponding to the bound constraints. If one of pwc or pwb
% is not passed to this routine, then no working set will be computed.

function [exitflag, x, y, z, wc, wb] = qppPostsolve(id, px, py, pz, pwc, pwb)

if (nargin < 4)
    error('qppPostsolve: Invalid number of input arguments!');
end

if (nargin < 6)
    pwb = [];
    if (nargin < 5)
        pwc = [];
    end
end


[exitflag, x, y, z, wc, wb] = qpp(11, id, px, py, pz, pwc, pwb, 0);

end
