function ccs_write_exe()
[n, m, ccs, i, v] = mat2ccs(tril([2,1,0,0,0;1,4,1,0,1;0,1,3,2,0;0,0,2,0,0;0,1,0,0,2]))
ccs_write(n, m, ccs, i, v)
ccs_write(n, m, ccs, i, v,'X.dat')

% example (matrix, rhs from hsl_ma77, sec5.1, first example). Write in ccsc format (format == 1)
[n, m, ccs, i, v] = mat2ccs([2, 3, 0, 0, 0; 3, 1, 4, 0, 6; 0, 4, 1, 5, 0; 0, 0, 5, 3, 0; 0, 6, 0, 0, 1])
ccs_write(n, m, ccs, i, v, [], 1)
ccs_write(n, m, ccs, i, v,'X.dat', 1)
