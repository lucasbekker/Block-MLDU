% ST3 loads the Stokes fluid flow problem system matrix X = [A, B''; B, 0]
% and the accompanying parameters. Run "ST_info" for more information. 

% Load the data.
load ST3.mat

% Print message.
fprintf('\nLoad: Stokes fluid flow problem system matrix X = [A, B''; B, 0]:\n');

% Extract parameters and matrices A, B. Include timing.
tic
n = nnz(diag(X)), N = size(X,1), m = N - n
A = X(1:n,1:n); B = X(n+1:N,1:n);
time = toc;

% Print message.
fprintf('Load: X is %d x %d; A is %d x %d; and B is %d x %d with m <=n\n',...
        (n + m),(n + m),n,n,m,n);

% Show timing results and clear.
time
clear time

% Provide a spy plot of X.
spy(X)
title({'Stokes fluid flow problem system matrix with spy()','  '});

% Print message to continue
fprintf('Press <ENTER> to continue ...\n');
pause

% Provide a spy plot of B.
spy(B)
title({'Its indicence matrix B', ' '});