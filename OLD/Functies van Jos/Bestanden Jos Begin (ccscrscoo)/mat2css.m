function [n, m, css, i, lu] = mat2css(A)
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v);
