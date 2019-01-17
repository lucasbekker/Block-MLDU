% truncate: leave sparse vectors 1:k
%      but: keep matrix dimensions
function [n, m, cxs, i, v] = ccs_truncate(n, m, cxs, i, v, k)
if k < length(cxs)
 cxs = cxs(1:k+1);
 i   = i(1:cxs(k+1)-1);
 v   = v(1:cxs(k+1)-1,:);
end
