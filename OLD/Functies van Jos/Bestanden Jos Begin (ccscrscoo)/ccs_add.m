function [n, m, ccs, i, v] = ccs_add(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, szr)
if nargin < 11 || isempty(szr), szr = 0; end;
[n, m, ccs, i, v] = ccs_fxy(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, @(x,y) x + y, szr);
