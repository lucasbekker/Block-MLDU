% ccs extensions are permitted, n := max(n, nk); m := max(m, k)
function [n, m, ccs, i, v] = ccs_set_sparse(n, m, ccs, i, v, k, nk, ik, vk, echo)
if (nargin < 10) || isempty(echo), echo = 0; end;
if echo > 0
 fprintf('Msg(ccs_set_sparse):\n');
end
if k > 0
 if k < length(ccs)
 % rewrite (insert, works also for ik, iv == [] (code uses length([]) == 0)
 % i(1:ccs(k)-1), ik, i(ccs(k+1):end)
  i = [i(1:ccs(k)-1); ik; i(ccs(k+1):end)];
  v = [v(1:ccs(k)-1,:); vk; v(ccs(k+1):end,:)];
  ccs = diff(ccs);
  ccs = [ccs(1:k-1); length(ik); ccs(k+1:end)];
  ccs = cumsum([1; ccs]);
 else
  % append empty columns
  ccs = [ccs; repmat(ccs(end), [k - length(ccs), 1])];
  % append column k
  [n, m, ccs, i, v] = ccs_stack_sparse(n, m, ccs, i, v, max(m, k), ik, vk);
 end
 n = max(n, nk);
else
 fprintf('Err(ccs_set_sparse): column index is negative: %d.\n',k);
end
