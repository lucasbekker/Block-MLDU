function css_ldu_exe(A, nb)
if nargin > 1
 css_ldu_example(A, nb);
else
fprintf('matrix 1 x 1 non-singular:\n');
css_ldu_example([2], '1')

fprintf('matrix 1 x 1 singular -- very special because its css representation is the empty (not zero!) css:\n');
fprintf(' factorization will indicate that this matrix is singular!\n');
css_ldu_example([0], '2')

fprintf('matrix 1 x 1 singular ... empty matrix:\n');
fprintf(' factorization is fine, factor L, D, U are all empty!\n');
css_ldu_example([], '3')

fprintf('matrix 2 x 2 non-singular is stored with css with 1 column ...:\n');
fprintf('css must extend during factorization:\n');
css_ldu_example([[2, -1]; [-1, 0]], '4')

fprintf('matrix 2 x 2 singular is stored with css with 1 column ...:\n');
css_ldu_example([[2, -1]; [0, 0]], '5')

fprintf('matrix 3 x 2 singular:\n');
css_ldu_example([[1, -1]; [1, 1]; [1, 2]], '6')

fprintf('matrix 2 x 3 singular:\n');
css_ldu_example([[1, -1, 1]; [-1, 0, 1]], '7')

fprintf('matrix 3 x 3 singular is stored with css with 1 column (first one)...:\n');
fprintf('css must extend during factorization:\n');
css_ldu_example([[1, -1, 1]; [-1, 0, 0]; [2, 0, 0]], '8')

fprintf('matrix 3 x 3 singular is stored with css with 1 column (middle one)...:\n');
fprintf(' singular with k == 0 valid columns in LU, so L, U and D are empty:\n');
css_ldu_example([[0, 0, 0]; [0, 1, -1]; [0, 1, 0]], '9')

fprintf('matrix full square 7 x 7 non-singular:\n');
css_ldu_example(magic(7), '10')

% example where L and U have integer entries, default point-wise factorization
% in addition, row/col 3 contains a zero diagonal entry, this entry is not stored in css (but the entries left and right of it are)...
% 1   0   1   0      1   0   1   0      1   0   1   0      1   0   1   0
% 2   1  -1   0 ->   2   1  -3   0 ->   2   1  -3   0  ->  2   1  -3   0
% 0  -1   0   1      0  -1   0   1      0  -1  -3   1      0  -1  -3   1
%-2   0   1   3     -2   0   3   3     -2   0   3   3     -2   0   3   4
%       A              GE row/col 1       GE row/col 2       GE row/col 3
css_ldu_example([1, 0, 1, 0; 2, 1, -1, 0; 0, -1, 0, 1; -2, 0, 1, 3], '11')

% example where L and U have integer entries, default point-wise factorization
% creates a temporary zero at the diagonal AND therefore a 0-width 4-th row/col lu vector ...
% 1   0   1   2      1   0   1   2      1   0   1   2      1   0   1   2
% 1   1   2   1 ->   1   1   1  -1 ->   1   1   1  -1  ->  1   1   1  -1      
% 0   1   2   1      0   1   2   1      0   1   1   2      0   1   1   2
% 3   4   4   2      3   4   1  -4      3   4  -3   0      3   4  -3   6
%      A              GE row/col 1       GE row/col 2       GE row/col 3
css_ldu_example([1, 0, 1, 2; 1, 1, 2, 1; 0, 1, 2, 1; 3, 4, 4, 2], '13')

% square 4 x 4 matrix with LU GE stops at update with row/col 3 (zero pivot)
% 1  0  1  2  -> 1  0  1  2  -> 1  0  1  2
% 1  1  2  1     1  1  1 -1     1  1  1 -1
% 0  1  1  1     0  1  1  1     0  1  0  2
% 3  4  4  2     3  4  1 -4     3  4 -3  0
%    A           GE row/col1   GE row/col2 --> (3,3) contains pivot 0 (stop)
css_ldu_example([1, 0, 1, 2; 1, 1, 2, 1; 0, 1, 1, 1; 3, 4, 4, 2], '14')

% rectangular 4 x 3 matrix with LU GE stops at update with row/col 3
% matrix is as in prev 4 x 4 example but without the last column
css_ldu_example([1, 0, 1; 1, 1, 2; 0, 1, 1; 3, 4, 4], '15')

% rectangular 3 x 4 matrix
% 1  0  1  2  -> 1  0  1  2  -> 1  0  1  2
% 1  1  2  1     1  1  1 -1     1  1  1 -1
% 0  1  2  1     0  1  2  1     0  1  1  2
css_ldu_example([1, 0, 1, 2; 1, 1, 2, 1; 0, 1, 2, 1], '16')

end

function css_ldu_example(A, nb)
epsilon = 2^-42;
fprintf('Msg(css_lcu_exe): STARTING test %s:\n',nb);
full_A = full(A)
fprintf('Msg(css_lcu_exe): Above matrix converted to css:\n',nb);
[n, m, css, i, lu] = mat2css(A); n, m , css_tp = css', i_tp = i', lu_tp = lu'
[n, m, css, i, lu] = css_ldu(n, m, css, i, lu); n, m, css_tp = css', i_tp = i', lu_tp = lu'
tic; LU = css2mat(n, m, css, i, lu); fprintf('Msg(css_lcu_exe): Factorization time: %f\n',toc); full_LU = full(LU)
% if amount of stored row/columns is less than min(n,m) then automatically, there is a zero pivot (since k < n, a_{nn} is not stored ...)
% whether this is ok (non-square A to start with) or not (singular A, where non-singular expected) is up to the user to decide
k = length(css) - 1;
n,m,k
if k < min(n,m) fprintf('Msg(css_lcu_exe): %sA. PASSED: lcu stopped on zero pivot at entry (%d, %d) -- matrix not square or singular!\n',nb, k+1, k+1);
elseif n == m
 fprintf('Test ALL (only if amount of stored columns in css equals %d = k == min(n, m)):\n', k);
 tic; L = tril(LU); D = diag(diag(LU,0)); U = triu(LU); fprintf('Msg(css_lcu_exe): Extraction time: %f\n',toc);
 tic; norm(A - L*(D\U), 'fro'); fprintf('Msg(css_lcu_exe): Norm calculation time: %f\n', toc);
 if norm(A - L*(D\U), 'fro') > epsilon , fprintf('Msg(css_lcu_exe): %sB. FAILED\n',nb); else fprintf('Msg(css_lcu_exe): %sB. PASSED\n',nb); end
end
fprintf('Test VER for A(:,1:%d) amount of columns:\n', k);
LUV = LU(:,1:k); tic; L = mtril(LUV); D = mdiag(LUV); U = mtriu(LUV); fprintf('Msg(css_lcu_exe): Extraction time: %f\n',toc);
full_L = full(L), full_D = full(D), full_U = full(U)
tic; norm(A(:,1:k) - L*(D\U), 'fro'); fprintf('Msg(css_lcu_exe): Norm calculation time: %f\n', toc);
if norm(A(:,1:k) - L*(D\U), 'fro') > epsilon , fprintf('Msg(css_lcu_exe): %sC. FAILED\n',nb); else fprintf('Msg(css_lcu_exe): %sC. PASSED\n',nb); end
fprintf('Test HOR for A(1:%d,:) amount of rows:\n', k);
LUH = LU(1:k,:); L = mtril(LUH); D = mdiag(LUH); U = mtriu(LUH);
if norm(A(1:k,:) - L*(D\U), 'fro') > epsilon , fprintf('Msg(css_lcu_exe): %sD. FAILED\n',nb); else fprintf('Msg(css_lcu_exe): %sD. PASSED\n',nb); end
