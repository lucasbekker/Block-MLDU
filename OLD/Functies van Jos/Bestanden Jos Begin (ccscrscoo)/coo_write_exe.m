function coo_write_exe()
% examples
fprintf('Msg(coo_write_exe): Example 1: No check\n');
A = [2,1,0,0,0;1,4,1,0,1;0,1,3,2,0;0,0,2,0,0;0,1,0,0,2];
[n, m, i, j, v] = mat2coo(A);
coo_write(n, m, i, j, v);
coo_write(n, m, i, j, v, 'A.dat');

% first example HSL2013 MA57
fprintf('Msg(coo_write_exe): Example 2: No check\n');
A = tril([2, 3, 0, 0, 0; 3, 0, 4, 0, 6; 0, 4, 1, 5, 0;0, 0, 5, 0, 0; 0, 6, 0, 0, 1])
[n, m, i, j, v] = mat2coo(A);
coo_write(n, m, i, j, v)
coo_write(n, m, i, j, v, 'A.dat');
