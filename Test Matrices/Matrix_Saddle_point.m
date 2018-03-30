function [ X ] = Matrix_Saddle_point(n,d)
% Relevant saddle point problem

% Input parsing
if nargin == 1
    d = 0.05;
end

% Parameters
N = n^2;
m = floor(N/3);

% Construct block matrices
A = Matrix_A(n);
B = triu(sprand(m,N,d),1) + [speye(m),sparse(m,N - m)];

% Construct final matrix
X = [A,B';B,sparse(m,m)];

end