% % l = sparse_find(k, i, r, l) let find the largest position p in 0:|i| such that i(p) < k
% %  return: 0 <= p <= |i| ==> i(1:p) < k and i(p+1:end) >= k
% % l, r are the initial positions used for a bisection search -- defaults are 0 and n+1
% % if there exists a q such that (given) i(q+1:end) > k then set 
%
function l = sparse_find(k, i, l, r)
if (nargin < 3) || isempty(l), l = 0; end; % assume sparse array begins at 0 with index 0 (never used)
if (nargin < 4) || isempty(r), r = length(i)+1; end; % assume sparse array ends at length(i)+1 with index k+1 (never used)
%fprintf('Attempt to find location of key %d in array(%d) ',k, length(i)); if size(i,2) > 1, i, else i_tp = i', end;
%fprintf('Using l = %d, r = %d\n',l, r);
% find the last index i <= k
while  (r > l + 1)
 m = floor((l + r)/2);
 if (i(m) < k), l = m; else r = m; end
end
