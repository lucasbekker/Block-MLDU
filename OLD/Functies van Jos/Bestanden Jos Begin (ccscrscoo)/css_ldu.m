% % factor A = (L+D) * D^{-1} * (U+D) returns LU where L = tril(LU), U = triu(LU) and D = diag(diag(LU)).
% % no pivotting, i.e., no permutation -- routine assumes A is positive definite, which implies all pivots
% % will be non-zero (can be shown). If a zero pivot is encountered in step k (which implies that A was not pos def)
% % then matrix part ccs(1:k) is returned -- that is the largest part of A for which
% % a non-singular factorization LD^-1U existed. Note: if A = [0] (1 x 1) then k = -1 on return!
% % 
% %
% Note, we calculate a Schur update for all 1x1 blocks b = 1,...,length(ccs)-1: A Schur update consists of:
%  1 Calculate D = ldu
%  2 Scale L and U with ldu
%  3 Calculate S = S - L d^-1 U
% For the last block b == length(s)-1 steps 1 ... 3 are not needed to obtain the factorization.
% We still execute 1 ... 3 and from 1 get information about whether the last diagonal block could be
% factored (and therefore guaranteed non-singular) -- step 2 and 3 are skipped.
% This permits us to return the largest ensured non-singular LDU factors.
% "Largest" in the sense that if a larger part of the input matrix would be factored,
% this factorization would be ensured to be singular.
% % 
% %
% % input n x m matrix A
% % succeeded (in factorizing the min(n,m) x min(n,m) part of A) if on return length(css) - 1 = min(n, m).
% %  failed if length(css) - 1 < min(n, m)
% %
% % (**) -- (see below): if m > n then css storage does not satisfy the sparse vector property [n, i, v] that n >= max(i)
% %                       for each column of the matrix, in fact columns v_i = [l_i, u_i] tend to be of length m > n
% %                       The write ccs_set_sparse(ccs matrix, sparse column of length m) would alter n -> m
% %                       which is not desirable. This is avoided by the line(s) (**) below
% %
% % remark: an (L+D)D^-1(U+D) factorization requires non-zero diagonal D to be stored -- this at its turn
% %          implies that all diagonal entries of the min(n,m) x min(n,m) part of A should be non-zero and had
% %          to be stored which at its turn implies that length(css) - 1 = min(n, m)
% %          upon a succesfull (partial) factorization.
% %
function [n, m, css, i, lu] = css_ldu(n, m, css, i, lu)
% extend css to ensure all row/col are stored, so access to row/col i(ik(d)) via css is ensures
% extended css is OUTPUT css, so is ok. Without this there should be extension
% inside the for d=1:length(ik) loop, which decreases change of parallellization of this loop
fprintf('*** LDU START:\n');
mx = max(n,m); % lu is (max(n,m) x 2)
mn = min(n,m); % css should have mn + 1 entries
epsilon = 2^-42;
k = 1;
% k < mn below is enough for LU but not enough for L * D^-1 * U because it is possible that the last entry is equal to 0
% in which case we would return a zero pivot (in D) and still have a failing factorization L * D^-1 * U
% note: for LU -- even if (only) U(end:end) = 0 -- one has L*U = A. This does not hold for L * D^-1 * U
% 
% so, for L * D^-1 * U we must process also the last entry (block) to determine whether it is singular
while (k < length(css)) % note: length(css) can increase inside the loop
 if (css(k+1) > css(k)) && (i(css(k)) == k) && (abs(lu(css(k),2)) > epsilon) % can Gauss Elim if diagonal entry exists and is non-zero
  fprintf('*** Work with diagonal block %d:\n', k);
  fprintf('*** Matrix (Schur complement) A^%d:\n', k-1); disp(full(css2mat(n, m, css, i, lu)))

  % factor d^{-1} = 1/d_kk: multiply l- or u-vector with it, here u-vector
  [nk, ik, luk] = ccs_get_sparse(n, m, css, i, lu, k);
  luk(:,2) = luk(:,2) / luk(1,2);
  % diagonal element no longer needed
  ik = ik(2:end);
  luk = luk(2:end, :);
  p = 1;
  fprintf('*** with row %d update cols/rows:\n',k); ldu_print_ik = disp(ik);
  ldu_print_luk = luk
  ldu_print_nk = nk


  %  for all entries in ik (d <= length(ik)) as long as row/col ik(d) is stored ...

  % With row/col k, eliminate row/col (ik(d))_{d = 2 ... min(n, m)} only since:
  % 0.  we assume that ik(1) =  k, the row to eliminate "with"
  % I.  only for ik(d) < n it is possible to eliminate an entry at row ik(d):
  %      for m > n, possibly ik(d) > n -- for which A_{ik(d),*} by definition == 0,
  %      (only (since m > n) A_{*,ik(d)} <> 0), so there is nothing to eliminate
  % II. only for ik(d) < m it is possible to eliminate an entry at row/col ik(d):
  %      if ik(d) > m then extension would add extra zero columns to the matrix ... non-sense
  while (p <= length(ik)) && (ik(p) <= mn)
   %fprintf('*** INSIDE loop for ik(p) = %d\n',ik(p));
   [nlu, ilu, vlu] = ccs_get_sparse(n, m, css, i, lu, ik(p));
   %fprintf('*** with   row/col %d with dofs > diagonal:\n',k); ik(p:end)'

   % choice right now for css is that diag(l) = diag(u), which requires double storage
   % and double calculation -- but makes the code more forward and straight.

   %% METHOD I: slower -- some folds could be skipped if luk(p,1) or luk(p,2) == 0
   %[nl, il, vl] = sparse_fxy(mx, ilu, vlu(:, 1), mx, ik(p:end), luk(p:end, 1), @(x,y) x - luk(p,2)*y);
   %[nu, iu, vu] = sparse_fxy(mx, ilu, vlu(:, 2), mx, ik(p:end), luk(p:end, 2), @(x,y) x - luk(p,1)*y);
   %[ns, is, lus] = sparse_join(nl, il, vl, nu, iu, vu);
   %[n, m, css, i, lu] = ccs_set_sparse(n, m, css, i, lu, ik(p), n, is, lus); %! use n i.o. ns because of (**)


   % METHOD II: faster -- if there is an optimal implementation behind (col-scale, only if scalefactor <> 0)
   lukd = - luk(p:end,:)*diag([luk(p,2), luk(p,1)]);
   ikd  = ik(p:end);
   [ns, is, lus] = sparse_fxy(mx, ilu, vlu, mx, ikd, lukd, @(x,y) x + y);
   [n, m, css, i, lu] = ccs_set_sparse(n, m, css, i, lu, ik(p), n, is, lus); %! use n i.o. ns because of (**)
   
   p = p + 1;
  end
 else
  fprintf('*** zero pivot at diagonal %d\n',k);
 
  % The returned part contains all entries of (k-1) x m UNION n x (k-1)
  % Gaussian Elimination created Schur complement: Factorization is
  % L * D^{-1} * U where D contains the pivots, so the last step (#k) which
  % created a zero (singular) block at the end of D should be excluded.
  % User can decide to use one or both of the factorizations
  %(k-1) x m matrix LU(1:k-1,:) or
  % n x (k-1) matrix LU(:,1:k-1) 
  css = css(1:k);
  i =   i(1:css(k)-1);
  lu = lu(1:css(k)-1,:);

  % force a stop
  k = k + 1;
 end

 k = k + 1;
end
fprintf('*** LDU END:\n');
