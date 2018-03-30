function [ E_A, E_X ] = Test_Function_1( n, t )
% This test function constructs the A matrix and the associated saddle
% point problem matrix. It subsequently calculates the block MLDU
% factorization and the Frobenius norm of the error.

% Parameters
N = n^2;

% Generate matrices
A = Matrix_A(n);
X = Matrix_Saddle_Point(n);

% Generate blocks and permutation
s = 2*ones(N/2,1);
S = [2*ones(1,m),ones(1,(N - m))];                                         % m 2x2 blocks and (N - m) 1x1 blocks
P = [reshape([1:m;(N + 1):(N + m)],[1,2*m]),(m + 1):N];

% Apply permutation
Y = X(P,P);

if ~t
    % Calculate Block MLDU factorizations (no timings)
    [L_A,D_A,U_A] = MLDU(A,s);
    [L_Y,D_Y,U_Y] = MLDU(Y,S);
else
    % Calculate Block MLDU factorizations (with timings)
    tic
    [L_A,D_A,U_A] = MLDU(A,s);
    toc
    tic
    [L_Y,D_Y,U_Y] = MLDU(Y,S);
    toc
end

% Calculate error norms
E_A = norm(A - (L_A + D_A)*(D_A\(D_A + U_A)),'fro');
E_X = norm(Y - (L_Y + D_Y)*(D_Y\(D_Y + U_Y)),'fro');

end