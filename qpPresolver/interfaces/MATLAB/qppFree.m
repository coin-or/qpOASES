%% qppFree
%  Destroy qpPresolver entity by deallocating memory
%---------------------------------------------------------------------------------
%
% Call qppFree( id ) to free the presolver entity with the given id.

function [] = qppFree(id)

if (nargin < 1)
    error('qppFree: Invalid number of input arguments!');
end

qpp(20, id);

end
