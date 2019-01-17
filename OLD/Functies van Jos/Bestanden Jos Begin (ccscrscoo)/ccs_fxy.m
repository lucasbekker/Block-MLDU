% apply a sparse fxy to the columns of ccs1 and ccs2, need not have same amount of columns
% too few columns will be solved by substitution zero columns for ccs1 or ccs2
function [n, m, ccs, i, v] = ccs_fxy(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, f, szr)
if nargin < 12 || isempty(szr), szr = 0; end;

[n, m, ccs, i, v] = ccs_zeros(max(n1,n2), max(m1,m2));

% add columns mutually stored in A and B
lc1 = length(ccs1)-1; lc2 = length(ccs2)-1;
mn = min(lc1, lc2)
for k=1:mn
 [nk, ik, vk] = sparse_fxy(n1, i1([ccs1(k):(ccs1(k+1)-1)]),v1([ccs1(k):(ccs1(k+1)-1)]), n2, i2([ccs2(k):(ccs2(k+1)-1)]),v2([ccs2(k):(ccs2(k+1)-1)]), f, szr);
 [n, m, ccs, i, v] = ccs_stack_sparse(n, m, ccs, i, v, nk, ik, vk);
end
if lc1 > mn
 empty_val2 = sparse(0, size(v2,2)); % empty col of ccs2 should have empty value of size (0 x .)
 for k=mn+1:lc1
  [nk, ik, vk] = sparse_fxy(n1, i1([ccs1(k):(ccs1(k+1)-1)]),v1([ccs1(k):(ccs1(k+1)-1)]), n2, [], empty_val2, f, szr);
 [n, m, ccs, i, v] = ccs_stack_sparse(n, m, ccs, i, v, nk, ik, vk);
 end
elseif lc2 > mn
 empty_val1 = sparse(0, size(v1,2)); % empty col of ccs2 should have empty value of size (0 x .)
 for k=mn+1:lc2
  [nk, ik, vk] = sparse_fxy(n1, [], empty_val1, n2, i2([ccs2(k):(ccs2(k+1)-1)]),v2([ccs2(k):(ccs2(k+1)-1)]), f, szr);
  [n, m, ccs, i, v] = ccs_stack_sparse(n, m, ccs, i, v, nk, ik, vk);
 end
end
