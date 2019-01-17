% % Much less efficient than crs_T_inv_crs_B due to k^2 scattered adaptations rather than about k times k grouped adaptations
% % Assume that matrix T
% %  - k x k
% %  - triangular 
% % Cases: T = L lower triangular and B is k x n, k >= 1, T=L stored ccs, B stored crs, output crs X = inv(L)*B
% %        T = U upper triangular and B is n x k, k >= 1, T=U stored crs, B stored ccs, output ccs X = B*inv(U)
% %               (this situation occurs when [LU] are stored css and one needs to compute B*inv(U) := L (d+u)^-1
% %                with B == L and U = d+u)
% % One implemenations does 2 different multiplications which can not be told apart from the inputs
% % nT, mT, ccsT, iT, vT, mBT, nBT, ccsBT, iBT, vBT The only way to tell them apart is to tell this routine
% % which version it calculates (is needed for echo only), so best is to switch versions on echo = 1, 2.
%
% % Assumes: T = L low and B crs: |ccsT| = |ccsBT|
% %          T = U upp and B ccs: |crsT| = |ccsB|
% % Assumes: i(c[cr]s(k)) = k is a stored non-zero diagonal entry
% % Assumes: T(i(c[cr]s(k))) <> 0 is in addition non-zero
%
function [mBT, nBT, ccsBT, iBT, vBT] = cxs_T_inv_cxs_B(nT, mT, ccsT, iT, vT, mBT, nBT, ccsBT, iBT, vBT, store_zero_results, echo)
if (nargin < 11) || isempty(store_zero_results), store_zero_results = 0; end;
if (nargin < 12) || isempty(echo), echo = 0; end;

if     (echo == 1)
 fprintf('Msg(cxs_T_inv_cxs_B): Input T UPPER is %d x %d:\n', nT, mT); disp(full(ccs2mat(nT, mT, ccsT, iT, vT))');
 fprintf('Msg(cxs_T_inv_cxs_B): Input B       is %d x %d:\n', mBT, nBT); disp(full(ccs2mat(mBT, nBT, ccsBT, iBT, vBT)));
elseif (echo == 2)
 fprintf('Msg(cxs_T_inv_cxs_B): Input T LOWER is %d x %d:\n', nT, mT); disp(full(ccs2mat(nT, mT, ccsT, iT, vT)));
 fprintf('Msg(cxs_T_inv_cxs_B): Input B       is %d x %d:\n', nBT, mBT); disp(full(ccs2mat(mBT, nBT, ccsBT, iBT, vBT))');
end

% We explicity proceed over all entries in T in order i(1:length(i)), so ccs of T is implicitly used
if (nT > 0) && abs(nT - mT) == 0 && abs((length(ccsT)-1) - (length(ccsBT)-1)) == 0
 for k=1:length(ccsT)-1

  % zero-rows/columns need not stored in ccs/crs format for B
  % scaling such rows/columns leads to still the zero rows/columns
  % and scaled-substracting them from the other columns does alter the other columns
  if ccsBT(k+1) > ccsBT(k)

   % scale row/col k of B (can be done in situ in Fortran/C)
   [nBk, iBk, vBk] = ccs_get_sparse(mBT, nBT, ccsBT, iBT, vBT, k, echo);
   fprintf('B row|col %d /:= %f\n', iT(ccsT(k)), vT(ccsT(k)));
   vBk = vBk / vT(ccsT(k));
   [mBT, nBT, ccsBT, iBT, vBT] =  ccs_set_sparse(mBT, nBT, ccsBT, iBT, vBT, k, nBk, iBk, vBk);

   % subtract row/col B_k from row/col r of B (can be done in situ in Fortran/C)
   for c=ccsT(k)+1:ccsT(k+1)-1
    r  = iT(c);
    [nBr, iBr, vBr] = ccs_get_sparse(mBT, nBT, ccsBT, iBT, vBT, r, echo);
    fprintf('B row %d -:= row|col %d * %f\n', c, r, vT(c));
    [nBr, iBr, vBr] =  sparse_fxy(nBr, iBr, vBr, nBk, iBk, vBk, @(x,y) x - vT(c)*y, store_zero_results);
    [mBT, nBT, ccsBT, iBT, vBT] =  ccs_set_sparse(mBT, nBT, ccsBT, iBT, vBT, r, nBr, iBr, vBr);
   end
  end;


 end
else
 if nT <= 0, fprintf('Msg(ccs_T_inv_crs_B): FAILED: nT = %d <= 0\n',nT); end;
 if abs(nT - mT) > 0, fprintf('Msg(ccs_T_inv_crs_B): FAILED: nT = %d <> %d = mT\n',nT,mT); end;
 if abs(length(ccsT) - (length(ccsBT))) > 0, fprintf('Msg(ccs_T_inv_crs_B): FAILED: |ccsT|-1 = %d <> %d = |ccsBT|-1\n',length(ccsT)-1,length(ccsBT)-1); end;
 abort;
end

