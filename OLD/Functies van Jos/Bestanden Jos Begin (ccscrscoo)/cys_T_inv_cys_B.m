% %
%
%  ONLY DIFFERENCE IN BEHAVIOUR WITH cxs version is that cys operates on the css of U' (and cxs on the css of L)
%  SEEMS THAT SAME INVARIANTS (MULTIPLICATIONS HOLD, see cxs_ and cys_ exe.m versions!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
%
%
%
% % Assume that matrix T
% %  - n x n
% %  - triangular 
% % Cases:
% %        T = U upper triangular and B is n x k, k >= 1, U stored crs, B stored crs, output crs X = inv(U)*B
% %        T = L lower triangular and B is k x n, k >= 1, L stored ccs, B stored ccs, output ccs X = B*inv(L)
% % One implemenations does 2 different multiplications -- if input matrix dimensions are swapped
% % (input dimensions need to be swap for T = L case because calls to routine cxs_get_col(), i suspect)
%
% % Assumes: i(c[cr]s(k)) = k is a stored non-zero diagonal entry
% % Assumes: T(i(c[cr]s(k))) <> 0 is in addition non-zero
%
function [nB, mB, crsB, iB, vB] = cys_T_inv_cys_B(nT, mT, crsT, iT, vT, nB, mB, crsB, iB, vB, store_zero_results, echo)
if (nargin < 11) || isempty(store_zero_results), store_zero_results = 0; end;
if (nargin < 12) || isempty(echo), echo = 0; end;
if abs(nT - mT) > 0, fprintf('FAILED: nT = %d <> %d = mT\n',nT,mT); end;
if abs(length(crsT) - (length(crsB))) > 0, fprintf('Msg(ccs_T_inv_crs_B): FAILED: |ccsT|-1 = %d <> %d = |crsB|-1\n',length(ccsT)-1,length(crsB)-1); end;

if     (echo == 1)
 fprintf('Msg(cxs_T_inv_cxs_B): Input T UPPER is %d x %d:\n', nT, mT); disp(full(ccs2mat(nT, mT, ccsT, iT, vT))');
 fprintf('Msg(cxs_T_inv_cxs_B): Input B       is %d x %d:\n', mBT, nBT); disp(full(ccs2mat(mBT, nBT, ccsBT, iBT, vBT)));
elseif (echo == 2)
 fprintf('Msg(cxs_T_inv_cxs_B): Input T LOWER is %d x %d:\n', nT, mT); disp(full(ccs2mat(nT, mT, ccsT, iT, vT)));
 fprintf('Msg(cxs_T_inv_cxs_B): Input B       is %d x %d:\n', nBT, mBT); disp(full(ccs2mat(mBT, nBT, ccsBT, iBT, vBT))');
end

% We explicity proceed over all entries in U in order i(1:length(i)), so crs of U is implicitly used
if abs(nT - mT) == 0 && abs((length(crsT)-1) - (length(crsB)-1)) == 0
 for k=length(crsT)-1:-1:1
  [nBk, iBk, vBk] = ccs_get_sparse(nB, mB, crsB, iB, vB, k);
  for c=crsT(k)+1:crsT(k+1)-1
   l  = iT(c);
   fprintf('B row %d -:= row %d * %f\n', k, l, vT(c));
   [nBl, iBl, vBl] = ccs_get_sparse(nB, mB, crsB, iB, vB, l);
   [nBk, iBk, vBk] = sparse_fxy(nBk, iBk, vBk, nBl, iBl, vBl, @(x,y) x - vT(c)*y, store_zero_results);
  end
  fprintf('B row %d /:= %f\n', iT(crsT(k)),  vT(crsT(k)));
  vBk = vBk / vT(crsT(k));
  [nB, mB, crsB, iB, vB] = ccs_set_sparse(nB, mB, crsB, iB, vB, k, nBk, iBk, vBk);
 end
else
 crsB = [1]; iB = []; vB= []; 
end

