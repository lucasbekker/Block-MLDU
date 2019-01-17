function [ E ] = Test_Function_1( n, MLDU_Function )
% This test function constructs the A matrix and subsequently calculates
% the block MLDU factorization and the Frobenius norm of the error.

% Parameters
N = n^2;

% Generate matrix
A = Matrix_A(n);

% Generate blocks
s = 2*ones(N/2,1);

% Calculate Block MLDU factorization
[L,D,U] = MLDU_Function(A,s);

% Calculate error norm
E = norm(A - (L + D)*(D\(D + U)),'fro');

end