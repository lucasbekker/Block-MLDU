function [ E ] = Test_Function_2( n, MLDU_Function )
% This test function constructs the saddle point problem matrix. It 
% subsequently calculates the block MLDU factorization and the Frobenius
% norm of the error.

% Parameters
N = n^2;
m = floor(N/3);

% Generate matrix
X = Matrix_Saddle_Point(n);

% Generate blocks and permutation
S = [2*ones(1,m),ones(1,(N - m))];                                         % m 2x2 blocks and (N - m) 1x1 blocks
P = [reshape([1:m;(N + 1):(N + m)],[1,2*m]),(m + 1):N];

% Apply permutation
Y = X(P,P);

% Calculate Block MLDU factorization
[L,D,U] = MLDU_Function(Y,S);

% Calculate error norms
E = norm(Y - (L + D)*(D\(D + U)),'fro');

end