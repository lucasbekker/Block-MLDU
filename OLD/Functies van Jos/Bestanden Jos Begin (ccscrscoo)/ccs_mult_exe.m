function ccs_mult_exe()
% Assume I(:,k) contains N(:,k) non-zero row numbers for column k (and perhaps a tail of fill-up zeros)
% Then copying all these row number vectors into one column vector for ccs can be done in matlab:
clear, clc
I = [1, 2, 5;3, 4, 0; 0, 7, 0]
i = I(find(I))
if max(abs(i - [1; 3; 2; 4; 7; 5])) > 0, fprintf('Msg(ccs_mult_exe): 1. FAILED\n'); else fprintf('Msg(ccs_mult_exe): 1. PASSED\n'); end;

% % incompatible matrices -> [] product. CAN NOT USE ccs_mult_exe (since A*B does not exist)
A = [1]
B = magic(3)
full(A*B)
[n1, m1, ccs1, i1, v1] = mat2ccs(A);
[n2, m2, ccs2, i2, v2] = mat2ccs(B);
[n, m, ccs, i, v] = ccs_mult(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, [], 0);
C = ccs2mat(n, m, ccs, i, v); full(C)
if norm([]-C, 'fro') > 0, fprintf('Msg(ccs_mult_exe): 2A. FAILED\n'); else fprintf('Msg(ccs_mult_exe): 2A. PASSED\n');  end;

% matlab 2014b, octave 3.8: sparse(2,0)*sparse(0,3) -> sparse(2,3) matrix
ccs_mult_example(sparse(2,0),sparse(0,3), '2');

% matlab 2014b, octave 3.8: sparse(0,3)*sparse(3,2) -> sparse(0,2) matrix
% i.e., matrix with dimension zero * non-zero matrix is matrix with dimension zero
A = sparse(0,3);
B = magic(3); B = B(:,[1,2]);
ccs_mult_example(A, B, '3');

% matlab 2014b, octave 3.8: sparse(4,3)*sparse(3,0) -> sparse(4,0) matrix
% i.e., matrix with dimension zero * non-zero matrix is matrix with dimension zero
A = magic(4); A = A(:,1:3);
B = sparse(3,0);
ccs_mult_example(A, B, '4');

% nil potent matrix product with itself is zero matrix, stores zero entry
% when sparse_add does not compress out zero entries (results of "a + (-a)")
A = [0, 1; 0, 0];
B = A;
ccs_mult_example(A, B,'5');

% matrices A and B with zero columns at end which are not stored
A = magic(4); A(1:4,3:4) = 0
B = magic(4); B(1:4,4) = 0
ccs_mult_example(A, B,'5X');

ccs_mult_example([ 1,  0; 0,  1],[-2, -1; 0, -3], '6');

A = toeplitz([1;0;0;2;3;0;0;4;0]); B = A;
ccs_mult_example(A, B,'7');

A = toeplitz([1;0;0;2;0;3;0;4;0],[1;5;0;6;0;0;7;8;0]);
B = A([4,6,2,3,8,7,1,9,5],[4,6,2,3,8,7,1,9,5]); B = fliplr(B);
ccs_mult_example(A, B,'8');

% diagonal band matrix A and its inverse B, A*B = I should create 5+4+...+1 = 15 zero entries
A = toeplitz([1,-1,0,0,0,0],[1,0,0,0,0,0]);
B = toeplitz([1,1,1,1,1,1],[1,0,0,0,0,0]);
ccs_mult_example(A, B,'9');

% nil potent matrix product with itself is zero matrix, stores zero entry
% when sparse_add does not compress out zero entries (results of "a + (-a)")
A = toeplitz([0,1,0,0,0,0,0],zeros(1,7))
B = A
full_A_power_5 = full(A*A*A*A*A)
[n, m, ccs, i, v] = mat2ccs(A);
[n2, m2, ccs2, i2, v2] = mat2ccs(B);
[n, m, ccs, i, v] = ccs_mult(n, m, ccs, i, v, n2, m2, ccs2, i2, v2); % A*A
fprintf('Amount of stored      entries of product C = A^2: %d.\n',ccs(end)-1);
fprintf('Amount of stored zero entries of product C = A^2: %d.\n',sum(v==0));
[n, m, ccs, i, v] = ccs_mult(n, m, ccs, i, v, n2, m2, ccs2, i2, v2); % A^3
fprintf('Amount of stored      entries of product C = A^3: %d.\n',ccs(end)-1);
fprintf('Amount of stored zero entries of product C = A^3: %d.\n',sum(v==0));
[n, m, ccs, i, v] = ccs_mult(n, m, ccs, i, v, n2, m2, ccs2, i2, v2); % A^4
fprintf('Amount of stored      entries of product C = A^4: %d.\n',ccs(end)-1);
fprintf('Amount of stored zero entries of product C = A^4: %d.\n',sum(v==0));
[n, m, ccs, i, v] = ccs_mult(n, m, ccs, i, v, n2, m2, ccs2, i2, v2); % A^5
fprintf('Amount of stored      entries of product C = A^5: %d.\n',ccs(end)-1);
fprintf('Amount of stored zero entries of product C = A^5: %d.\n',sum(v==0));
C = ccs2mat(n, m, ccs, i, v); full_A_power_5 = full(C)
if norm(A*A*A*A*A-C, 'fro') > 0, fprintf('Msg(ccs_mult_exe): 10. FAILED\n'); else fprintf('Msg(ccs_mult_exe): 10. PASSED\n');  end;

function ccs_mult_example(A, B, nb)
if abs(size(A,2) - size(B,1)) > 0
 fprintf('Err(ccs_mult_example): incompatible dimensions!\n');
 abort
end
full_A = full(A), full_B = full(B), full_A_times_B = full(A*B)
[n1, m1, ccs1, i1, v1] = mat2ccs(A);
[n2, m2, ccs2, i2, v2] = mat2ccs(B);
[n, m, ccs, i, v] = ccs_mult(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2);
fprintf('Amount of stored      entries of product C = A*B: %d.\n',ccs(end)-1);
fprintf('Amount of stored zero entries of product C = A*B: %d.\n',sum(v==0));
C = ccs2mat(n, m, ccs, i, v); full_A_times_B = full(C)
A_times_B_exact = full(A*B)
if norm(A_times_B_exact-C, 'fro') > 0, fprintf('Msg(ccs_mult_exe): %sA. FAILED\n',nb); else fprintf('Msg(ccs_mult_exe): %sA. PASSED\n',nb); end;
%% test adjoint operation
%[m1, n1, ccs1T, i1T, v1T] = mat2ccs(A');
%[m2, n2, ccs2T, i2T, v2T] = mat2ccs(B');
%[m, n, ccsT, iT, vT] = ccs_mult(m1, n1, ccs1T, i1T, v1T, m2, n2, ccs2T, i2T, v2T);
%fprintf('Amount of stored      entries of product C = B*A: %d.\n',ccs(end)-1);
%fprintf('Amount of stored zero entries of product C = B*A: %d.\n',sum(v==0));
%C = ccs2mat(m, n, ccsT, iT, vT)'; full_B_times_A = full(C)
%B_times_A_exact = full(B*A)
%if norm(B_times_A_exact-C, 'fro') > 0, fprintf('Msg(ccs_mult_exe): %sB. FAILED\n',nb); else fprintf('Msg(ccs_mult_exe): %sB. PASSED\n',nb); end;
