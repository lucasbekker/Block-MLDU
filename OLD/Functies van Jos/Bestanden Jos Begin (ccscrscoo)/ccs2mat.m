function A = ccs2mat(n, m, ccs, i, v, echo)
if (nargin < 6) || isempty(echo), echo = 0; end;
if echo
 fprintf('Msg(ccs2mat): Input matrix (ccs format):\n');
 ccs_write(j, m, ccs, i, v, [], echo);
end
% due too circum stances, css stored matrices store can zero values outside the n x m matrix
% filter out these entries first, the non-zero entries must fit into the n x m matrix
% matlab sparse checks whether all (i,j) are inside (1...n,1...m):
[n, m, i, j, v] = ccs2coo(n, m, ccs, i, v);
s = (abs(v)>0);
A = sparse(i(s),j(s),v(s),n,m); % full(U)
