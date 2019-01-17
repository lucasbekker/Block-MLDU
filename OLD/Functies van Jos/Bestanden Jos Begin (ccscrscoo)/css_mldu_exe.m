% Add example which fails in second pivot block s=[1,2,1] and see whether size of the largest minor returned == 2


function css_mldu_exe(A, s, nb)
if nargin > 1
 css_mldu_example(A, s, nb)
else
fprintf('matrix 1 x 1 non-singular:\n');
css_mldu_example([2], [], '1')

fprintf('matrix 1 x 1 singular -- very special because its css representation is the empty (not zero!) css:\n');
fprintf(' factorization will indicate that this matrix is singular!\n');
fprintf('matrix has only non-singular 0 x 0 diagonal block:\n');
css_mldu_example([0], [], '2', 0)

fprintf('matrix 1 x 1 empty matrix -- SHOULD NOT BE FLAGGED SINGULAR, i.e., k=n=m=0 should be upon return:\n');
fprintf(' factorization is fine, factor L, D, U are all empty!\n');
css_mldu_example([], [], '3')

fprintf('matrix 2 x 2 non-singular is stored with css with 1 column ...:\n');
fprintf('css must extend during factorization:\n');
css_mldu_example([[2, -1]; [-1, 0]], [], '4')

fprintf('matrix 2 x 2 singular is stored with css with 1 column ...:\n');
fprintf('matrix has only non-singular 1 x 1 diagonal block:\n');
css_mldu_example([[2, -1]; [0, 0]], [], '5', 1)

fprintf('matrix 3 x 2 singular:\n');
css_mldu_example([[1, -1]; [1, 1]; [1, 2]], [], '6')

fprintf('matrix 2 x 3 singular:\n');
css_mldu_example([[1, -1, 1]; [-1, 0, 1]], [], '7')

fprintf('matrix 3 x 3 singular is stored with css with 1 column ...:\n');
fprintf('css must extend during factorization:\n');
fprintf('matrix has only non-singular 2 x 2 diagonal block:\n');
css_mldu_example([[1, -1, 1]; [-1, 0, 0]; [2, 0, 0]], [], '8', 2)

fprintf('matrix 3 x 3 singular is stored with css with 1 column (middle one)...:\n');
fprintf(' singular with k == 0 valid columns in LU, so L, U and D are empty:\n');
fprintf('matrix has only non-singular 0 x 0 diagonal block:\n');
css_mldu_example([[0, 0, 0]; [0, 1, -1]; [0, 1, 0]], [], '9', 0)

fprintf('matrix full square 7 x 7 non-singular:\n');
css_mldu_example(magic(7), [], '10')

fprintf('[1, 2] blocks on a 2 x 3 matrix, 2-nd diagonal block partially outside the matrix:\n');
css_mldu_example([2, -1, 0;-1, 2, 3], [1, 2], '10INF');


% example where L and U have integer entries, default point-wise factorization
% in addition, row/col 3 contains a zero diagonal entry, this entry is not stored in css (but the entries left and right of it are)...
% 1   0   1   0      1   0   1   0      1   0   1   0      1   0   1   0
% 2   1  -1   0 ->   2   1  -3   0 ->   2   1  -3   0  ->  2   1  -3   0
% 0  -1   0   1      0  -1   0   1      0  -1  -3   1      0  -1  -3   1
%-2   0   1   3     -2   0   3   3     -2   0   3   3     -2   0   3   4
%       A              GE row/col 1       GE row/col 2       GE row/col 3
css_mldu_example([1, 0, 1, 0; 2, 1, -1, 0; 0, -1, 0, 1; -2, 0, 1, 3], [], '11')

% example where L and U have integer entries, default point-wise factorization
% creates a temporary zero at the diagonal AND therefore a 0-width 4-th row/col lu vector ...
% 1   0   1   2      1   0   1   2      1   0   1   2      1   0   1   2
% 1   1   2   1 ->   1   1   1  -1 ->   1   1   1  -1  ->  1   1   1  -1      
% 0   1   2   1      0   1   2   1      0   1   1   2      0   1   1   2
% 3   4   4   2      3   4   1  -4      3   4  -3   0      3   4  -3   6
%      A              GE row/col 1       GE row/col 2       GE row/col 3
css_mldu_example([1, 0, 1, 2; 1, 1, 2, 1; 0, 1, 2, 1; 3, 4, 4, 2], [], '13')

% square 4 x 4 matrix with LU GE stops at update with row/col 3 (zero pivot)
% 1  0  1  2  -> 1  0  1  2  -> 1  0  1  2
% 1  1  2  1     1  1  1 -1     1  1  1 -1
% 0  1  1  1     0  1  1  1     0  1  0  2
% 3  4  4  2     3  4  1 -4     3  4 -3  0
%    A           GE row/col1   GE row/col2 --> (3,3) contains pivot 0 (stop)
fprintf('matrix 4 x 4 with non-singular 2 x 2 diagonal block:\n');
css_mldu_example([1, 0, 1, 2; 1, 1, 2, 1; 0, 1, 1, 1; 3, 4, 4, 2], [], '14', 2)

% rectangular 4 x 3 matrix with LU GE stops at update with row/col 3
% matrix is as in prev 4 x 4 example but without the last column
fprintf('matrix 4 x 3 with non-singular 2 x 2 diagonal block:\n');
css_mldu_example([1, 0, 1; 1, 1, 2; 0, 1, 1; 3, 4, 4], [], '15', 2)

% rectangular 3 x 4 matrix
% 1  0  1  2  -> 1  0  1  2  -> 1  0  1  2
% 1  1  2  1     1  1  1 -1     1  1  1 -1
% 0  1  2  1     0  1  2  1     0  1  1  2
css_mldu_example([1, 0, 1, 2; 1, 1, 2, 1; 0, 1, 2, 1], [], '16')

% example where L and U have integer entries, now super-node factorization. 22 is equal to 11
% s = [1, 1, 1, 1]
% 1   0   1   0  1   1   0   1   0  1   1   0   1   0   1  1   0   1   0
% 2   1  -1   0 ->   2   1  -3   0 ->   2   1  -3   0  ->  2   1  -3   0
% 0  -1   0   1      0  -1   0   1      0  -1  -3   1      0  -1  -3   1
%-2   0   1   3     -2   0   3   3     -2   0   3   3     -2   0   3   4
css_mldu_example([1, 0, 1, 0; 2, 1, -1, 0; 0, -1, 0, 1; -2, 0, 1, 3], [], '22')
% s = [2, 1, 1]
% 1   0   1   0  2   1   0   1   0  1   1   0   1   0
% 2   1  -1   0 ->   2   1  -1   0 ->   2   1  -1   0
% 0  -1  -3   1      0  -1  -3   1      0  -1  -3   1
%-2   0   3   3     -2   0   3   4     -2   0   3   4
%
% s = [1, 2, 1]
% 1   0   1   0  1   1   0   1   0  2   1   0   1   0
% 2   1  -1   0 ->   2   1  -3   0 ->   2   1  -3   0
% 0  -1  -3   1      0  -1   0   1      0  -1   0   1
%-2   0   3   3     -2   0   3   3     -2   0   3   4
css_mldu_example([1, 0, 1, 0; 2, 1, -1, 0; 0, -1, 0, 1; -2, 0, 1, 3], [2, 1, 1], '23')
css_mldu_example([1, 0, 1, 0; 2, 1, -1, 0; 0, -1, 0, 1; -2, 0, 1, 3], [2, 2], '24')
css_mldu_example([1, 0, 1, 0; 2, 1, -1, 0; 0, -1, 0, 1; -2, 0, 1, 3], [1, 2, 1], '25')

fprintf('matrix full square 7 x 7 non-singular:\n');
css_mldu_example(magic(7), [2, 1, 3, 1], '101')
fprintf('matrix full square 11 x 11 non-singular, but second pivot block starts with 10^-14, so rank 2:\n');
fprintf('matlab rank command returns 11.\n');
css_mldu_example(magic(11), [2, 3, 1, 3, 2], '102', 2)

fprintf('matrix sparse square 9 x 9 non-singular block toeplitz: scalar (1x1 block) factorization:\n');
css_mldu_example(toeplitz_2d_2wn(3, [-1, -1, 4, -1, -1]), [], '201') % point-wise factorization
fprintf('matrix sparse square 9 x 9 non-singular block toeplitz: block factorization:\n');
css_mldu_example(toeplitz_2d_2wn(3, [-1, -1, 4, -1, -1]), [3,3,3], '202') % block-factorization
fprintf('matrix sparse square 36 x 36 non-singular block toeplitz: scalar (1x1 block) factorization:\n');
css_mldu_example(toeplitz_2d_2wn(6, [-1, -1, 4, -1, -1]), [], '203') % point-wise factorization
fprintf('matrix sparse square 36 x 36 non-singular block toeplitz: block factorization:\n');
css_mldu_example(toeplitz_2d_2wn(6, [-1, -1, 4, -1, -1]), [6,6,6,6,6,6], '204') % block-factorization
fprintf('matrix sparse square 36 x 36 non-singular block toeplitz: [3,3,3] block factorization -- will be finished with 1x1 blocks:\n');
css_mldu_example(toeplitz_2d_2wn(6, [-1, -1, 4, -1, -1]), [3,3,3], '205') % block-factorization -- mistake in block size input

fprintf('matrix sparse square non-symmetric and positive definite because eig((A + A'')/2) in (0,infty), scalar factorization:\n');
css_mldu_example(toeplitz_2d_2wn(6,[-1, -1, 4, -1, -1/2]), [], '210') % point-wise factorization
fprintf('matrix sparse square non-symmetric and positive definite because eig((A + A'')/2) in (0,infty), block factorization:\n');
css_mldu_example(toeplitz_2d_2wn(6,[-1, -1, 4, -1, -1/2]), [1, 2, 3, 2, 2], '211') % block factorization

fprintf('matrix X sparse square 9 x 9 non-singular constraint with matrix A 6 x 6 positive definite: scalar (1x1 block) factorization:\n');
A = toeplitz([2,-1,0,0,0,0]), A = sparse(A);
B = speye(3,6); B(1,3) = -1; B(2,3) = -1; B(3,6) = -1; full(B)
X = [A, B'; B, sparse(3,3)]; full(X)
css_mldu_example(X, [], '300') % point-wise factorization, can fail ...
fprintf('matrix X sparse square 9 x 9 non-singular constraint with matrix A 6 x 6 positive definite: block [2, 2, 2, 1, 1, 1] factorization:\n');
p = [1, 6, 2, 7, 3, 9, 4, 5, 6]; Y=X(p,p)
css_mldu_example(X, [2, 2, 2], '301') % block factorization of permuted matrix, must work


end

function css_mldu_example(A, s, nb, k_nsng)
if nargin < 4 || isempty(k_nsng), k_nsng = min(size(A)); end; % default assumption is that s x s diagonal is non-singular (s = min(n,m))
epsilon = 2^-36;
fprintf('Msg(css_mldu_exe): EXAMPLE %s:\n',nb);
full_A = full(A)
fprintf('Msg(css_mldu_exe): Above matrix converted to css:\n',nb);
fprintf('Msg(css_mldu_exe): Input block size sequence s:\n',nb); disp(s);
%s_input = s;
[n, m, css, i, lu] = mat2css(A); n, m , css_tp = css', i_tp = i', lu_tp = lu'
[n, m, css, i, lu, cssX, iX, luX] = css_mldu(n, m, css, i, lu, s); n, m, css_tp = css', i_tp = i', lu_tp = lu'
%[n, m, css, i, lu, cssX, iX, luX, s] = css_mldu(n, m, css, i, lu, s); n, m, css_tp = css', i_tp = i', lu_tp = lu'
%fprintf('Msg(css_mldu_exe): Input block size sequence s:\n',nb); disp(s_input);
fprintf('Msg(css_mldu_exe): Output block size sequence s:\n',nb); disp(s);
fprintf('VERDICT:\n');
tic; LU = css2mat(n, m, css, i, lu); fprintf('Msg(css_mldu_exe): Factorization time: %f\n',toc); full_LU = full(LU)
tic; LUX = css2mat(n, m, cssX, iX, luX); fprintf('Msg(css_mldu_exe): Factorization time: %f\n',toc); full_LU = full(LU)
% if amount of stored row/columns is less than min(n,m) then automatically, there is a zero pivot (since k < n, a_{nn} is not stored ...)
% whether this is ok (non-square A to start with) or not (singular A, where non-singular expected) is up to the user to decide
k = length(css) - 1;
n,m,k
mn = min(n,m);
if k < mn
 % not maximal rank case
 if k == k_nsng
  fprintf('Msg(css_mldu_exe): %s S. PASSED: mldu: Largest non-singular block of input %d x %d matrix is %d x %d\n',nb, n, m, k, k);
  fprintf('Msg(css_mldu_exe): %s S. PASSED: mldu: MATRIX IS SINGULAR AND NOT OF MAX RANK\n',nb);
 else
  fprintf('Msg(css_mldu_exe): %s S. FAILED: mldu: Largest non-singular block of input %d x %d matrix is %d x %d\n',nb, n, m, k, k);
  fprintf('Msg(css_mldu_exe): %s S. FAILED: mldu: MATRIX IS SINGULAR AND NOT OF MAX RANK\n',nb);
  fprintf('Msg(css_mldu_exe): %s S. FAILED: mldu: Flagged as failure because predicted largest non-singular block was %d x %d block\n',nb, k_nsng, k_nsng)
 end
else
 % maximal rank case
 if n ~= m
  fprintf('Msg(css_mldu_exe): %s M. PASSED: mldu: The non-square %d x %d matrix is of maximal rank\n',nb, n, m);
  fprintf('Msg(css_mldu_exe): %s M. PASSED: mldu: MATRIX NON-SQUARE AND MAXIMAL RANK\n',nb);
 else
  fprintf('Msg(css_mldu_exe): %s N. PASSED: mldu: The %d x %d matrix is square and non-singular\n',nb, n, m);
  fprintf('Msg(css_mldu_exe): %s N. PASSED: mldu: MATRIX SQUARE AND NON-SINGULAR!\n',nb);
 end
end
fprintf('Test BLOCK-triangular LDU for A(:,1:%d) -- VERTICAL:\n', k);
LUV = LU(:,1:k); tic; L = mtril(LUV,0,s); D = mdiag(LUV,0,s); U = mtriu(LUV,0,s); fprintf('Msg(css_mldu_exe): Extraction time: %f\n',toc);
if prod(size(A)) < 100
 full_L = full(L), full_D = full(D), full_U = full(U)
end
tic; norm_difference_BLOCK_triangular = norm(A(:,1:k) - L*(D\U), 'fro'), fprintf('Msg(css_mldu_exe): Norm calculation time: %f\n', toc);
if norm(A(:,1:k) - L*(D\U), 'fro') > epsilon , fprintf('Msg(css_mldu_exe): %s BLOCK-triangular V. FAILED\n',nb); else fprintf('Msg(css_mldu_exe): %s BLOCK-triangular V. PASSED\n',nb); end
fprintf('Test BLOCK-triangular LDU for A(1:%d,:) -- HORIZONTAL:\n', k);
LUH = LU(1:k,:); L = mtril(LUH,0,s); D = mdiag(LUH,0,s); U = mtriu(LUH,0,s);
if prod(size(A)) < 100
 full_L = full(L), full_D = full(D), full_U = full(U)
end
if norm(A(1:k,:) - L*(D\U), 'fro') > epsilon , fprintf('Msg(css_mldu_exe): %s BLOCK-triangular H. FAILED\n',nb); else fprintf('Msg(css_mldu_exe): %s BLOCK-triangular H. PASSED\n',nb); end
fprintf('Test TRIANGULAR LDU for A(:,1:%d) -- VERTICAL:\n', k);
LUV = LUX(:,1:k); tic; L = mtril(LUV,0); D = mdiag(LUV,0); U = mtriu(LUV,0); fprintf('Msg(css_mldu_exe): Extraction time: %f\n',toc);
if prod(size(A)) < 100
 full_L = full(L), full_D = full(D), full_U = full(U)
end
tic; norm_difference_TRIANGULAR = norm(A(:,1:k) - L*(D\U), 'fro'), fprintf('Msg(css_mldu_exe): Norm calculation time: %f\n', toc);
if norm(A(:,1:k) - L*(D\U), 'fro') > epsilon , fprintf('Msg(css_mldu_exe): %s TRIANGULAR V. FAILED\n',nb); else fprintf('Msg(css_mldu_exe): %s TRIANGULAR V. PASSED\n',nb); end
fprintf('Test TRIANGULAR LDU for A(1:%d,:) -- HORIZONTAL:\n', k);
LUH = LUX(1:k,:); L = mtril(LUH,0); D = mdiag(LUH,0); U = mtriu(LUH,0);
if prod(size(A)) < 100
 full_L = full(L), full_D = full(D), full_U = full(U)
end
if norm(A(1:k,:) - L*(D\U), 'fro') > epsilon , fprintf('Msg(css_mldu_exe): %s TRIANGULAR H. FAILED\n',nb); else fprintf('Msg(css_mldu_exe): %s TRIANGULAR H. PASSED\n',nb); end
