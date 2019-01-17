% sparse matrix multiplication for two sparse matrices A, B in ccs format. Let C = A*B.
% calculates all columns C_l = sum_k A_k b_{kl} where A_k is column k of A and b_{kl} is non-zero
% entry k,l of B. C_l += A_k b_{kl} is calculated with sparse_fxy().
%
% efficient implementation:
%  - amount of index searches is the amount done by sparse_fxy()
%  - ccs data for C is obtained by simple stack, no search or sort needed
%
%
function [n, m, ccs, i, v] = ccs_mult(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, szr, echo)
if nargin < 11 || isempty(szr), szr = 0; end;
if nargin < 12 || isempty(echo), echo = 1; end;

[n, m, ccs, i, v] = ccs_zeros(0, 0);

if abs(m1 - n2) > 0
 if echo > 0
  fprintf('Err(ccs_mult): Matrix dimensions incompatible --> return [] matrix.\n');
 end
else
 % calculate column l of A*B for l=1:size(B,2). Each column is a linear combination of columns of A
 n = n1; m = m2;
 for l=1:length(ccs2)-1 % calculate column l of A*B which is [A*B]_l = sum_{i=1} -- columns of B without entries are not stored and lead to zero columns anyway, so don't need to be considered
  [nl, il, vl] = sparse_zeros(0);
  for k=ccs2(l):ccs2(l+1)-1
   beta = v2(k);
   r = i2(k);

   % columns r with r >= |ccs1| are not stored and are zero so don't contribute
   if r < length(ccs1)
    [n1s, i1s, v1s] = ccs_get_sparse(n1, m1, ccs1, i1, v1, r)
    [nl, il, vl] = sparse_fxy(nl, il, vl, n1s, i1s, v1s, @(x,y)x + beta*y, szr);
   end

  end
  [n, m, ccs, i, v] = ccs_stack_sparse(n, m, ccs, i, v, nl, il, vl);
 end
end
