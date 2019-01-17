function [n, m, ccs, i, v] = ccs_stack_sparse(n, m, ccs, i, v, nk, ik, vk, echo)
if (nargin < 9) || isempty(echo); echo = 0; end;
if echo > 0, fprintf('Msg(ccs_stack_sparse):\n'); end;
n = max(n, nk); % sparse array [nk, ik, vk] represents sequence in R^{nk} and we assume that sum(abs(v(i>n))) == 0
m = max(m, length(ccs));
i = [i; ik];
v = [v; vk];
ccs = [ccs; ccs(end) + length(ik)]; % don't use size(i,1) (i can be 1 x 0)
