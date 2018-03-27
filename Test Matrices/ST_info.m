clear
fprintf('Information on stokes system matrices X = [A, B''; B, 0]:\n')
for i=1:5, filename = strcat('ST',strcat(sprintf('%1d',i),'.mat')); load(filename); fprintf('size of X of stokes system %s is %d x %d\n', filename, size(X,1), size(X,2)); end
fprintf('Extract A and B with:\n')
fprintf('n = nnz(diag(X)), N = size(X,1), m = N - n\n')
fprintf('A = X(1:n,1:n); B = X(n+1:N,1:n);\n')
