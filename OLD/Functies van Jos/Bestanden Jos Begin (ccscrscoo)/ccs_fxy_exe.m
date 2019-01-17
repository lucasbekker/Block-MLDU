function ccs_fxy_exe()
clear, clc

% % matrices A and B with zero columns at end which are not stored
A = magic(4); A(1:4,3:4) = 0
B = magic(4); B(1:4,4) = 0
ccs_fxy_example(A, B, @(x,y) x + y, '1A');

% matrix addition of sparse and (fully stored) sparse matrix
A = toeplitz([-1;0;3;0;-5]);
B = magic(5);
ccs_fxy_example(A, B, @(x,y) x + y, '1');

% matrix addition of sparse A and -A leads to zero matrix
A = toeplitz([-1;0;3;0;-5]);
B = -A;
ccs_fxy_example(A, B, @(x,y) x + y, '2');

% matrix with mostly zero upper triangular part crs of that part only stores rows 1 & 2 whereas ccs stores columns 1:6 ...
A = tril(magic(7),-1), A(2,2) = 1
B = tril(magic(7),-1); B(2,2) = 1; B = B'
ccs_fxy_example(A, B, @(x,y) x + y, '3');

function ccs_fxy_example(A, B, f, nb)
fullA = full(A), fullB = full(B), full_A_plus_B_matlab_plus = full(A+B)
[n1, m1, ccs1, i1, v1] = mat2ccs(A);
[n2, m2, ccs2, i2, v2] = mat2ccs(B);
[n, m, ccs, i, v] = ccs_fxy(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, f);
fprintf('Amount of stored      entries of product C = A*B: %d.\n',ccs(end)-1);
fprintf('Amount of stored zero entries of product C = A*B: %d.\n',sum(v==0));
C = ccs2mat(n, m, ccs, i, v);
full_A_plus_B_fold_calculated = full(C)
if norm(C-(A+B), 'fro') > 0 || abs(n - n1) > 0 || abs(m - m2) > 0, fprintf('Msg(ccs_fxy_exe): %s. FAILED\n',nb); else fprintf('Msg(ccs_fxy_exe): %s. PASSED\n',nb); end;
