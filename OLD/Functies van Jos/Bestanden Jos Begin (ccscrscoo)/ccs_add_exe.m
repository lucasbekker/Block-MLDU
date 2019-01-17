function ccs_add_exe()
% Assume I(:,k) contains N(:,k) non-zero row numbers for column k (and perhaps a tail of fill-up zeros)
% Then copying all these row number vectors into one column vector for ccs can be done in matlab:
clear, clc

% % matrices A and B with zero columns at end which are not stored
A = magic(4); A(1:4,3:4) = 0
B = magic(4); B(1:4,4) = 0
ccs_add_example(A, B, '1');

% matrix addition of sparse and (fully stored) sparse matrix
A = toeplitz([-1;0;3;0;-5]);
B = magic(5);
ccs_add_example(A, B, '2');

% matrix addition of sparse A and -A leads to zero matrix
A = toeplitz([-1;0;3;0;-5]);
B = -A;
ccs_add_example(A, B, '3');

% matrix with mostly zero upper triangular part crs of that part only stores rows 1 & 2 whereas ccs stores columns 1:6 ...
A = tril(magic(7),-1), A(2,2) = 1
B = tril(magic(7),-1); B(2,2) = 1; B = B'
ccs_add_example(A, B, '4');

function ccs_add_example(A, B, nb)
fullA = full(A), fullB = full(B), full_A_plus_B = full(A+B)
[n1, m1, ccs1, i1, v1] = mat2ccs(A);
[n2, m2, ccs2, i2, v2] = mat2ccs(B);
[n, m, ccs, i, v] = ccs_add(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2);
fprintf('Amount of stored      entries of product C = A*B: %d.\n',ccs(end)-1);
fprintf('Amount of stored zero entries of product C = A*B: %d.\n',sum(v==0));
C = ccs2mat(n, m, ccs, i, v);
full_A_plus_B = full(C)
if norm(C-(A+B), 'fro') > 0 || abs(n - n1) > 0 || abs(m - m2) > 0, fprintf('Msg(ccs_add_exe): %s. FAILED\n',nb); else fprintf('Msg(ccs_add_exe): %s. PASSED\n',nb); end;
