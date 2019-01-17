function [n, m, ccs, i, v] = coo2ccs(n, m, i, j, v)
% permute [i, j, v] such that columns j are in ascending order
[j, p] = sort(j); % STATEMENT NOT NEEDED IF ONE ASSUMES coo TO BE j-sorted (as is demanded in readme.txt)
% REPLACE SORT BY OPTIONAL CHECK (DEFAULT CHECK: ON) IF SORT IS REMOVED
i = i(p);
v = v(p);
if length(j) > 0
 ccs = zeros(j(end),1);
 for c=1:length(j)
  ccs(j(c)) = ccs(j(c)) + 1;
 end
else
 ccs = [];
end 
ccs = cumsum([1;ccs]);
