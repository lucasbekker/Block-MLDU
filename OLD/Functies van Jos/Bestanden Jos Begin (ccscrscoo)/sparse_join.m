% % sparse_join joins two sparse columns (column matrices) into
% %                   one sparse column matrices
%
% % TYPICAL USE: (1) join two columns v and w into sparse U := [v, w]
% %              after which, one can implement
% %              (2) v+w, v+a*w on full arrays via U(1)+U(2), U(1)+a*U(2)
%
function [n, i, w] = sparse_join(n1, i1, v1, n2, i2, v2)
sparse_invariants(n1, i1, v1);
sparse_invariants(n2, i2, v2);


[n, i, w] = sparse_zeros(max(n1, n2));

n1 = length(i1); n2 = length(i2);
d1 = max(1,size(v1,2)); d2 = max(1,size(v2,2)); % join the [] vector into a ``n x 1'' vector (instead of into a ``n x 0'' vector (without min(1,.) construction)

l = 1; r = 1;
while  (l <= length(i1)) && (r <= length(i2))
 d = i1(l) - i2(r);
 if     d == 0
  [n, i, w] = sparse_stack(n, i, w, i1(l), [v1(l,:), v2(r,:)]);
  l = l+1; r = r+1; 
 elseif d > 0
  [n, i, w] = sparse_stack(n, i, w, i2(r), [zeros(1,d1), v2(r,:)]);
  r = r+1;
 else
  [n, i, w] = sparse_stack(n, i, w, i1(l), [v1(l,:), zeros(1,d2)]);
  l = l+1;
 end
end
% faster for non-fold case
if     l <= n1
 [n, i, w] = sparse_stack(n, i, w, i1(l:n1), [v1(l:n1,:), zeros(n1-l+1,d2)]);
elseif r <= n2
 [n, i, w] = sparse_stack(n, i, w, i2(r:n2), [zeros(n2-r+1,d1), v2(r:n2,:)]);
end
