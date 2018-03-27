clear;clc

load('S3D7u.mat')

S = [2*ones(1,m),ones(1,n-m)]; % m 2x2 blokken en N-m 1x1 blokken op de diagonal

P  = micro_blocks( n , m );
%P = [reshape([1:m;n+1:n+m],[1,2*m]),m+1:n];

Y = X(P,P);

[ L,D,U] = MLDU(Y,S);

norm(Y-(L+D)*(D\(D+U)),'fro')