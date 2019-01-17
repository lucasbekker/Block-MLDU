function A = coo2mat(n, m, i, j, v)
A = sparse(i, j, v, n, m);
