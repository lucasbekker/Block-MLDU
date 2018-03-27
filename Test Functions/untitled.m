clear, clc

n = 10; N = n^2;

A = kron(eye(n),gallery('tridiag',n,-1,2,-1))+kron(gallery('tridiag',n,-1,2,-1),eye(n));

[ L,D,U] = MLDU_test5(A,2*ones(N/2,1));
%[ L,D,U] = MLDU_3(A,2*ones(N/2,1));

norm(A-(L+D)*(D\(D+U)),'fro')

 

% % Relevant saddle point probleem: A sym pos def zoals bv hierboven en B zoals beneden
% 
% d = 0.05;
% 
% m = floor(N/3);
% 
% B = triu(sprand(m,N,d),1) + [speye(m),sparse(m,N-m)];
% 
% X = [A,B';B,sparse(m,m)];
% 
% S = [2*ones(1,m),ones(1,N-m)]; % m 2x2 blokken en N-m 1x1 blokken op de diagonal
% 
% P = [reshape([1:m;N+1:N+m],[1,2*m]),m+1:N];
% 
% Y = X(P,P);
% 
% [ L,D,U,S3,S4] = MLDU_test6(Y,S);
% 
% norm(Y-(L+D)*(D\(D+U)),'fro')