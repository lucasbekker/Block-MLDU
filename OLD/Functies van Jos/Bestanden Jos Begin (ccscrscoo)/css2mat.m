% % function css2mat returns either (1) entire matrix (nargout == 1) exor
% %  (2) its strictly lower triangular part LA and upper triangular part UA
% % ASSUMPTION: diag(L) = diag(U)
%
function [LU, U] = css2mat(n, m, css, i, lu, echo)
if nargin < 6 || isempty(echo), echo = 0; end;
% we use i = [] (i.e., i = 0 x 0 -- but we should have used i = 0 x 1, the empty vector),
% but since v = 0 x 1 abs(v) > 0 & i ... fails (because objects of 0 x 0 and 0 x 1 are compared for the empty matrix)
if length(i) > 0


 %[n, m, i, j, v] = ccs2coo(n, m, css, i, lu(:,1));
 %% if sparse(.,n,m) does not work then non-zero entries were outside n x m block which is forbidden
 %s = (abs(v)>0) & (i ~= j); % strip on abs(v) > 0 to delete zero entries L_{ij} in L with i > n, elsewise next sparse() fails for n < m
 %LU = sparse(i(s),j(s),v(s),n,m); % full(LU)
 LU = ccs2mat(n, m, css, i, lu(:,1), echo);
 LU = LU - spdiags(diag(LU), 0, n, m);
 %LU(speye(n,m)==1) = [];

 %[n, m, i, j, v] = ccs2coo(n, m, css, i, lu(:,2));
 %s = (abs(v)>0); % U should only have zero entries outside n x m block, so sparse(.,n,m) below should work
 %U = sparse(i(s),j(s),v(s),m,n)'; % full(U)
 U = ccs2mat(m, n, css, i, lu(:,2), echo)';

else
 LU = sparse(n,m);
 U = sparse(n,m);
end
if nargout < 2, LU = LU + U; end
