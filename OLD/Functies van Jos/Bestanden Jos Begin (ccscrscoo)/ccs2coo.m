function [n, m, i, j, v] = ccs2coo(n, m, ccs, i, v)
j = []; for k=1:length(ccs)-1, j = [j; repmat(k, [ccs(k+1)-ccs(k), 1])]; end;
