% simpler example:
clear, clc
N = 4;
A = spdiags([-ones(N,1),2*ones(N,1),-ones(N,1)],[-1,0,1],N,N);
[ L,D,U] = MLDU(A,2*ones(N/2,1));
norm(A-(L+D)*(D\(D+U)),'fro')


% Relevant saddle point probleem: A sym pos def zoals bv hierboven en B zoals beneden

d = 0.3;

m = floor(N/3);

B = [1,1,0,0];
%B = triu(sprand(m,N,d),1) + [speye(m),sparse(m,N-m)];

X = [A,B';B,sparse(m,m)];

S = [2*ones(m,1);ones(N-m,1)]; % m 2x2 blokken en N-m 1x1 blokken op de diagonal
%S = ones(5,1)
P = [reshape([1:m;N+1:N+m],[1,2*m]),m+1:N];

Y = X(P,P);
full(Y)

[ L,D,U] = MLDU(Y,S);
% U = U(:,(1:end-1));
% L = L((1:end-1),:);
% full(L)
% full(D)
% full(U)

norm(Y-(L+D)*(D\(D+U)),'fro')