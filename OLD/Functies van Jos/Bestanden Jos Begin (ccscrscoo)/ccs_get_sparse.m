% implementation:
%  simulates that all m columns of n x m ccs are stored
%  does not alter stored data (only read access)
%
% does not check whether returned max(i) > n ...
% will be the case, but such that v_{k:k>n} = 0 for all k
%
function [no, io, vo] = ccs_get_sparse(n, m, ccs, i, v, k, echo)
if (nargin < 7) || isempty(echo), echo = 0; end;
if echo > 0
 fprintf('Msg(ccs_get_sparse): Input matrix (ccs format):\n');
 ccs_write(n, m, ccs, i, v, [], [], echo);
end
if (k > 0) && (k <= m)
 start = ccs(end); stop = ccs(end) - 1;
 if k < length(ccs)
  start = ccs(k); stop = ccs(k+1)-1;
 end
 io = i(start:stop);
 vo = v(start:stop, :);
 no = n;
else
 fprintf('Err(ccs_get_sparse): column index out of range: %d not inside [1, %d].\n',k, m); abort;
end
