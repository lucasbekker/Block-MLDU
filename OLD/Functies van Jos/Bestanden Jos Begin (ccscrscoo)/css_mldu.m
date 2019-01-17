% % factor a n x m positive definite complex valued matrix A in css format without pivotting return both
% %
% % input: n x m complex valued matrix A in css format
% % input: n integer sequence s: sequence of pivot block sizes to be used for factorization
% % output:
% %
% % A = (Lb + Db)Db^-1(Db + Ub) with
% %   - Db       diagonal block-triangular
% %   - Lb strictly lower block-triangular
% %   - Ub strictly upper block-triangular
% %   - stored together in a single matrix [Lb,Db,Ub] in css format: [n, m, ccs, i, lu] (overwrites the input css A matrix)
% % A = (L + D)D^-1(U + D) with
% %   - D       diagonal triangular
% %   - L strictly lower triangular
% %   - stored together in a single matrix [L,D,U] in css format: [n, m, ccsLDU, iLDU, luLDU]
% %
% % Let ccsX be one of ccs, ccsLDU: The matrix A was of full rank if
% % |cssX| - 1 = min(n, m), and of course non-singular if additionally n = m.
% % If |ccsX| - 1 = k < min(n, m) then (due to no pivotting) only the first k rows/columns of A
% % were linearly independent and Lb*Db^{-1}*Ub resp. L*D^{-1}*U are equal to the matrix
% % which contains the first k rows and columns of A (remnant zero-filled).
% %
% % Multiplications with Db^-1 are involved in the calculation of the Schurcomplements (for Lb and Ub).
% % Because the end user (caller) also has to perform multiplications with Db^{-1} it is sensible
% % to return Db^{-1} related information (such as a factorization) to the caller of this routine.
% % Or to fully avoid such a return by returning
% %
% % [(Lb+Db)Db^-1, Db+Ub], where (Lb+Db)Db^-1 has implied non-saved unit-block-diagonal
% %
% % which is conform standard LU factorizations, but implies L <> U^T
% % for symmetric matrices A, undesirable.
% %
% % This leave us with approaches:
% % 1. - return BLOCK STRICTLY TRIANGULAR AND DIAGONAL parts [(Lb + Db), Db, (Db + Ub)] in a single matrix
% %    - next user has to calculate x -> Db^{-1} x
% %    - because this routine needs something similar to X -> Db^{-1} X to calculate
% %      the Schur update S = S - (Lb+Db)*Db^{-1}(Db+Ub) this approach is not optimal
% %      (user would need to for instance factor Db again)
% %      This routine uses ldu to factor Db which is ensured to be ok only for positive definite matrices ... (**)
% %      OBSERVATION: CALCULATING Db^{-1} WITH MATLAB GIVES MUCH LARGER ERRORS:
% %         ||A - (Lb + Db)Db^{-1}(Db + Ub)||_Fro = +- 100 ||A - L*D^{-1}*U||_Fro (below)!!
% %
% % 2. - return SCALAR STRICTLY TRIANGULAR and DIAGONAL parts [(L + D), D, (D + U)] in a single matrix
% %    - this scalar factorization (non-block one) can be obtained from [(Lb + Db), Db, (Db + Ub)] as follows:
% %    - non-pivot ldu factor Db = (l+d)d^-1(d+u) to calculate S = S - (Lb+Db)*Db^{-1}(Db+Ub) = S - [Db; Lb] Db^-1 [Db, Ub] (**)
% %    - where
% %         [Db; Lb] Db^-1 [Db, Ub] = 
% %         [(l+d)d^-1(d+u); L] (d+u)^-1 d (l+d)^-1 [(l+d)d^-1(d+u), U] = 
% %         [(l+d); L(d+u)^-1 d] d^-1 d (l+d)^-1 [(l+d)d^-1(d+u), U] = 
% %         [(l+d); L(d+u)^-1 d] d^-1  [(d+u), d (l+d)^-1U].
% %    - so the scalar factors (non-block) [(L + D), D, (D + U)] to be returned are
% %         [(l+d); L(d+u)^-1 d]: L
% %         [(d+u), d (l+d)^-1U]: U
% %         d^-1                : D^-1
% %      remark: To get a Choleski type factorization one would multipy both L and U with sqrt(d^-1)
% %              which is not possible if A is indefinite, so we don't attempt to do so.
% %
% %      (**) The non-pivot ldu factorization of Db must exist (for X = [A, B'; B, 0] it does exist if A is pos def and B = Schilders)
% %           and this does not hold for all non-singular Db (Db = [0, 1; 1, 0] does not
% %           have a ldu factorization (without pivotting)).
% %
% %      Technical remark: During 1 and 2 we calculate L(d+u)^-1 d.
% %      Remark: Making all columns of L have the same indices (and storing the extra zeros)
% %      is more expensive on storage than storing sparse L(d+u)^-1 d
% %      (because (d+u) is triangular), but perhaps (as in SUPERLU) it is cheaper in computer time
% %
% %      REMARK: For the Schilders B-orientation the factorization keeps the constraints intact because (Lb+Db)*Db^{-1}(Db+Ub) = L*D^{-1}*U
% %              Furthermore, a 2-column L block scaled by (u+d) keeps the triangular shape of B intact (is has a good change of doing so)
% %
% %
% %      REMARK: Below, we store U' in css format (same as U in crs format)
% %
% % examples:
% % (L+D)D^{-1}(D+U) of []    0x0 matrix result in kxk with k = -1 invertible part
% % (L+D)D^{-1}(D+U) of c=0   1x1 matrix result in kxk with k =  0 invertible part
% % (L+D)D^{-1}(D+U) of c<>0  1x1 matrix result in kxk with k =  1 invertible part: c = c * c^{-1} * c = l*d*u factorization
% %
function [n, m, css, i, lu, cssLdU, iLdU, luLdU] = css_mldu(n, m, css, i, lu, s)
if (nargin < 6) || isempty(s), s = []; end;
s = block_trim(cumsum([1, s]), min(n, m)) % also output s -- in case it needed to be adapted here, since user will use s for mtril, mtriu and mdiag
[n, m, cssLdU, iLdU, luLdU] = ccs_zeros(n, m);

% extend css to ensure all row/col are stored, so access to row/col i(ik(d)) via css is ensures
% extended css is OUTPUT css, so is ok. Without this there should be extension
% inside the for d=1:length(ik) loop, which decreases change of parallellization of this loop
mx = max(n,m);   % lu is (max(n,m) x 2)
mn = min(n,m);   % css should have mn + 1 entries
epsilon = 2^-42; % smaller numbers are zero
b = 1;           % start factorization with block 1


% Below, we calculate a Schur update for all blocks b = 1,..., length(s)-1: A Schur update consists of:
%  1 Calculate D = ldu
%  2 Scale L and U with ldu
%  3 Calculate S = S - L d^-1 U
% For the last block b == length(s)-1 steps 1 ... 3 are not needed to obtain the factorization.
% We still execute 1 ... 3 and from 1 get information about whether the last diagonal block could be
% factored (and therefore guaranteed non-singular) -- step 2 and 3 are skipped.
% This permits us to return the largest ensured non-singular LDU factors.
% "Largest" in the sense that if a larger part of the input matrix would be factored,
% this factorization would be ensured to be singular.
while (b < length(s)) % note: length(css) can increase inside the loop
 fprintf('Work with diagonal block %d:\n', b);
 fprintf('Entire matrix (includes Schur complement starting at entry (%d, %d)) A^%d (full format) is %d x % d:\n', s(b), s(b), b-1, n, m); disp(full(css2mat(n, m, css, i, lu)))


 % if the first row column of a micro block is empty or has the pivot not stored, or has zero pivot then we can not
 % determine the block-inverse with our ldu routine -- i.e., stop the factorization -- because ldu does not pivot
 % or matrix is singular
 k = s(b);
 if (k < length(css)) && (css(k+1) > css(k)) && (i(css(k)) == k) && (abs(lu(css(k),2)) > epsilon), cnt = true; else cnt = false; delta_nonzero_pivots = 0; end

 % cnt == true ==> D block exists and has a change to be non-singular (d_11 exists and d_11 <> 0). Collect the D block:
 if cnt == true
  % collect the D (css) and L,U (ccs) and U part related to the micro block
  [nD, mD, cssD, iD, vD] = ccs_zeros(s(b+1)-1, s(b+1)-1);
  fprintf('Diagonal block %d has order %d\n', b, nD-s(b)+1)
  [nLU, mLU, cssLU, iLU, vLU] = ccs_zeros(n, s(b+1)-s(b));

  while k < s(b+1)
   fprintf('Diagonal block %d contains row %d\n',b, k)
   % find in col/row k all entries < s(b+1) -- search with l = 0 and r - 1 = "max amount of possible entries s(b) <= i(p) < s(b+1)" which is s(b+1)-s(b)
   f = sparse_find(s(b+1), i(css(k):(css(k+1)-1)), 0, length(i(css(k):(css(k+1)-1)))+1);
   %fprintf('found location: %d\n',f);
   [nD, mD, cssD, iD, vD] = ccs_stack_sparse(nD, mD, cssD, iD, vD, nD, i(css(k):(css(k)+f-1)), lu(css(k):(css(k)+f-1),:));
   [nLU, mLU, cssLU, iLU, vLU] = ccs_stack_sparse(nLU, mLU, cssLU, iLU, vLU, n, i((css(k)+f):(css(k+1) - 1)), lu((css(k)+f):(css(k+1) - 1),:));
   k = k + 1;
  end % last value of k not used below

  % Shift D's dofs into 1, ..., s(b+1)-s(b):
  iD = iD - (s(b) - 1); nD = nD - (s(b) - 1); mD = mD - (s(b) - 1);

  fprintf('extracted D  block -- vars %d -- %d:\n', s(b), s(b+1)-1);
  full_extracted_D_block = full(css2mat(nD, mD, cssD, iD, vD))
  fprintf('extracted LU  block -- vars %d -- %d:\n', s(b), s(b+1)-1);
  fprintf('extracted LU  block can not be printed, is not a single matrix!\n');
  % instead of call to ldu, could use a recursive call to mldu (with s = [1,1,...,1])
  [nD, mD, cssDldu, iDldu, vDldu] = css_ldu(nD, mD, cssD, iD, vD);
  fprintf('Block D: D has been (L+D)D^-1(U+D) factored:\n');
  full_LCU_block = full(css2mat(nD, mD, cssDldu, iDldu, vDldu))


  delta_nonzero_pivots = length(cssDldu);
  if delta_nonzero_pivots < nD + 1
   cnt = false;
  end;
 end

  % If cnt == true then the pivot block D exists (!) and has been factorized into non-singular factors without pivoting.
  % Therefore invariants below after "cnt == true" are:
  % css(s(b)+1) > css(s(b)): first row column s(b) of block b contains an entry
  % i(css(s(b))) == s(b): diagonal entry D = d_kk exists (k := s(b)) (if it would not, D = 0 => D would be singular and factorization would have failed)
  % D = d_kk > epsilon (was tested for by css_ldu, because css_ldu does not pivot, |d_kk| <= epsilon, ==> ldu will stop with < maximal rank
  % So: mldu factorization of D exists



  if cnt == true



   % compute modified block parts [(l+d); L(d+u)^-1 d] d^-1  [(d+u), d (l+d)^-1U]
   %

  %fprintf('start to eliminate with row %d rows which follow below\n',k);

% % Below operation cxs_L_inv_cxs_B has to be used. This operation is less efficient than desired due to the
% % css storage of L. More efficient for cxs_L_inv_cxs_B would be a crs-like storage of the cols/rows of ldu:
% % +--------+
% % |\ ^^^^^^|
% % | \ ||||||
% % |<-\ |||||
% % |<--\ ||||
% % |<---\ |||
% % |<----\ ||
% % |<-----\ |
% % +--------+
% % Right now only less efficient css storage of ldu implemented (better put: ldu is output of css_ldu which output css format)


%
%    Note that the operations below can be done in-situ (i.e., directly applied to the LU block)
%


%
%    extract L block, is ccs stored
%
     [nL, mL, ccsL, iL, vL] = ccs_copy(n, mD, cssLU, iLU, vLU(:,1)); % iLU need not be shifted
     fprintf('extracted_L_block (full format) is %d x %d:\n', nL, mL); disp(full(ccs2mat(nL, mL, ccsL, iL, vL)))
     fprintf('extracted_L_block (css format):\n');
     fprintf('first line of output: n, m, #stored_cols+1, nnz:\n');
     ccs_write(nL, mL, ccsL, iL, vL);
%
%    extract U block, is crs stored, swap dimensions to have U' ccs stored (we only want ccs stored matrices)
%
     [mUT, nUT, ccsUT, iUT, vUT] = ccs_copy(m, nD, cssLU, iLU, vLU(:,2)); % iU need not be shifted
     fprintf('extracted_U_block (full format) is %d x %d:\n', nUT, mUT); disp(full(ccs2mat(mUT, nUT, ccsUT, iUT, vUT))')
     fprintf('extracted_U_T_block (ccs format) -- writes UT instead of U:\n');
     fprintf('%first line: n, m, #stored_rows+1, nnz:\n');
     ccs_write(mUT, nUT, ccsUT, iUT, vUT); % writes U^T

%
%    Check below can not be on cssLU because this stores both L and U
%
     if nL > 0 || mUT > 0
%
%    By performing matrix multiplications with flag "store non-zero results" == 1
%    we expect the sparsity patterns iUT and iL to be identical
%
      fprintf('calculate: L*inv(d+u):\n');
      [nL, mL, ccsL, iL, vL ] = cxs_T_inv_cxs_B(mD, nD, cssDldu, iDldu, vDldu(:,2), nL, mL, ccsL, iL, vL, 1, 1);
      fprintf('L*inv(d+u) (full format) is %d x %d:\n', nL, mL);
      full_L_invlpd = full(ccs2mat(nL, mL, ccsL, iL, vL))
      fprintf('amount of diagonal entries in d in ldu of D, and its diagonal entries, should all be non-zero:\n');
      nD, iDldu_tp = iDldu(cssDldu(1:end-1))', linsp_1_nD = 1:nD, vDldu_tp = vDldu(cssDldu(1:end-1))'
      fprintf('calculate: L*inv(d+u)*d:\n');
      [nL, mL, ccsL, iL, vL ] = ccs_scale_sparse(nL, mL, ccsL, iL, vL, nD, iDldu(cssDldu(1:end-1)), vDldu(cssDldu(1:end-1)));
      fprintf('L*inv(d+u)*d (full format) is %d x %d:\n', nL, mL);
      full_L_invlpd_d = full(ccs2mat(nL, mL, ccsL, iL, vL))

      fprintf('calculate: inv(l+d)*U:\n');
      [mUT, nUT, ccsUT, iUT, vUT] = cxs_T_inv_cxs_B(nD, mD, cssDldu, iDldu, vDldu(:,1), mUT, nUT, ccsUT, iUT, vUT, 1, 2);
      fprintf('inv(l+d)*U (full format) is %d x %d:\n', mUT, nUT);
      full_invlpd_U = full(ccs2mat(mUT, nUT, ccsUT, iUT, vUT)')
      fprintf('amount of diagonal entries in d in ldu of D, and its diagonal entries, should all be non-zero:\n');
      mD, iDldu_tp = iDldu(cssDldu(1:end-1))', linsp_1_mD = 1:mD, vDldu_tp = vDldu(cssDldu(1:end-1))'
      fprintf('calculate: d*inv(l+d)*U:\n');
      [mUT, nUT, ccsUT, iUT, vUT] = ccs_scale_sparse(mUT, nUT, ccsUT, iUT, vUT, mD, iDldu(cssDldu(1:end-1)), vDldu(cssDldu(1:end-1)));
      fprintf('d*inv(l+d)*U (full format) is %d x %d:\n', mUT, nUT);
      full_d_invlpd_U = full(ccs2mat(mUT, nUT, ccsUT, iUT, vUT)')
     else
      fprintf('L and U block are empty and are not scaled: L*(d+u)^-1*d respectively d*inv(l+d)*U\n');
     end

     if ~sequence_equal(ccsL, ccsUT, [], 1), fprintf('Err(css_mldu): FAILURE: Asymmetry in L and U sparsity patterns ccsL/ccsUT, should not happen ... (but not proven)\n'); end;
     if ~sequence_equal(iL, iUT, [], 1), fprintf('Err(css_mldu): FAILURE: Assymmetry in L and U sparsity patterns iL/iU, should not happen ... (but not proven)\n'); end;


  % Zip the scaled L and U blocks column/row by column/row into a css LU block
  for d=1:nD
   [nDx, iDx, vDx] = ccs_get_sparse(nD, nD, cssDldu, iDldu, vDldu, d);
   [nDlx, iDlx, vDlx] = ccs_get_sparse(nL, mL, ccsL, iL, vL, d);
   [nDux, iDux, vDux] = ccs_get_sparse(mUT, nUT, ccsUT, iUT, vUT, d);
   [n, m, cssLdU, iLdU, luLdU] = ccs_stack_sparse(n, m, cssLdU, iLdU, luLdU, n, [iDx + s(b) - 1; iDlx], [vDx; [vDlx, vDux]]);
  end

%
% Right now we need to compute S -= L U where L and U are as above (U' stored).
% Typically L, U^T  in \real^{n x d} whence L U = \sum_{l=1}^d L_l U_l
% where L_l and U_l^T are the column-l of L resp. U^T
%

  for d=1:nD
   fprintf('with row %d of Schur diagonal block %d -- %d prepare to update cols/rows:\n', s(b) + d - 1, s(b), s(b+1)-1);
   % factor d^{-1} = 1/d_kk: multiply l- or u-vector with it, here u-vector
   [nk, ik, lk] = ccs_get_sparse(nL, mL, ccsL, iL, vL, d);
   [mk, ik, uk] = ccs_get_sparse(mUT, nUT, ccsUT, iUT, vUT, d);
   % nk, mk not used
   diagonal_d_entry = vDldu(cssDldu(d))
   uk = uk / vDldu(cssDldu(d));

  % With [X:=L*(u+d)^-1*d,         Y:=d^-1 * d*(l+d)^-1*U] eliminate row/col (ik(p))_{p = 1 ... min(n, m)} only since:
  % I.  only for ik(p) <= n it is possible to eliminate an entry at row ik(p):
  %      for m > n, possibly ik(p) > n -- for which A_{ik(p),*} by definition == 0,
  %      (only (since m > n) A_{*,ik(p)} <> 0), so there is nothing to eliminate
  % II. only for ik(p) <= m it is possible to eliminate an entry at row/col ik(p):
  %      if ik(p) > m then extension would add extra zero columns to the matrix ... non-sense
  p = 1;
  fprintf('with row %d of Schur diagonal block %d -- %d prepare to update rows with dofs:\n', s(b) + d - 1, s(b), s(b+1)-1); disp(ik')
  while (p <= length(ik)) && (ik(p) <= mn)
   fprintf('with row %d of Schur diagonal block %d -- %d start update row:\n', s(b) + d - 1, s(b), s(b+1)-1); disp(ik(p))
   [nlu, ilu, vlu] = ccs_get_sparse(n, m, css, i, lu, ik(p));

   % choice right now for css is that diag(l) = diag(u), which requires double storage
   % and double calculation -- but makes the code more forward and straight.

   %% METHOD 0: slower -- some sparse_fxy could be skipped if luk(p,1) or luk(p,2) == 0
   %% WATCH OUT: fold does not store zero results default, so il <> iu possible below ==> sparse_join() is a must
   %[nl, il, vl] = sparse_fxy(n, ilu, vlu(:, 1), n, ik(p:end), lk(p:end), @(x,y) x - uk(p)*y);
   %[nu, iu, vu] = sparse_fxy(m, ilu, vlu(:, 2), m, ik(p:end), uk(p:end), @(x,y) x - lk(p)*y);
   %[ns, is, lus] = sparse_join(nl, il, vl, nu, iu, vu); % il <> iu possible ...
   %% END METHOD 0

   %% METHOD I: slower -- some sparse_fxy could be skipped if luk(p,1) or luk(p,2) == 0
   %% WATCH OUT: we here forse sparse_fxy to save zero results ==> il == iu and sparse_join() is not needed
   %[nl, il, vl] = sparse_fxy(n, ilu, vlu(:, 1), n, ik(p:end), lk(p:end), @(x,y) x - uk(p)*y, 1);
   %[nu, iu, vu] = sparse_fxy(m, ilu, vlu(:, 2), m, ik(p:end), uk(p:end), @(x,y) x - lk(p)*y, 1);
   %ns = mx; is = il; lus = [vl,vu];
   %% END METHOD I


   % METHOD II: faster -- if there is an optimal implementation behind (col-scale, only if scalefactor <> 0)
   % must use mx since both sequences of length n and m (in L and U parts) are combined into one sparse vector
   [ns, is, lus] = sparse_fxy(mx, ilu, vlu, mx, ik(p:end), [lk(p:end)*uk(p), uk(p:end)*lk(p)], @(x,y) x - y);
   % END METHOD II

   [n, m, css, i, lu] = ccs_set_sparse(n, m, css, i, lu, ik(p), n, is, lus); %! use n i.o. ns because of (**)
   
   p = p + 1;
  end
 end %loop over d=1:nD
 else
  % Typical case where ldu factorization fails is for A = [0] or A = [0, 1; 1, 0]
  fprintf('Scalar ldu factorization of %d x %d pivot block nb %d failed (could be singular, need not be)\n', s(b+1)-s(b), s(b+1)-s(b), b);
  fprintf('.... with zero pivot at position %d.\n', delta_nonzero_pivots + 1);

  % size of the largest non-singular minor is the sum of all sofar succesfully factored blocks
  % plus length(cssDldu)
  k = (s(b) - 1) + delta_nonzero_pivots;
 
  % The returned part contains all entries of (k-1) x m UNION n x (k-1)
  % Gaussian Elimination created Schur complement: Factorization is
  % L * D^{-1} * U where D contains the pivots, so the last step (#k) which
  % created a zero (singular) block at the end of D should be excluded.
  % User can decide to use one or both of the factorizations
  %(k-1) x m matrix LU(1:k-1,:) or
  % n x (k-1) matrix LU(:,1:k-1) 
  [n, m, css, i, lu] = ccs_truncate(n, m, css, i, lu, k);
  %css = css(1:k+1);
  %i =   i(1:css(k+1)-1);
  %lu = lu(1:css(k+1)-1,:);

  % force a stop
  b = length(s) + 1;
 end

 b = b + 1;
end
