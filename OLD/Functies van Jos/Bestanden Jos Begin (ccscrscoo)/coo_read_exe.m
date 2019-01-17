function coo_read_exe()
% example
A = [toeplitz([1, 0, -3, 0, 5, 0, 0, 7]),zeros(8,1)]
[n, m, i, j, v] = mat2coo(A);
coo_write(n, m, i, j, v, 'coo.data');
[n, m, i, j, v] = coo_read('coo.data');
B = coo2mat(n, m, i, j, v);
coo_read_example(A, B, '1');

A = stokes17();
[n, m, i, j, v] = mat2coo(A);
coo_write(n, m, i, j, v, 'coo.data');
[n, m, i, j, v] = coo_read('coo.data');
B = coo2mat(n, m, i, j, v);
coo_read_example(A, B, '2');

function coo_read_example(A, B, nb)
max_difference = max(max(abs(A - B))), if max_difference == 0, fprintf('Msg(coo_read_exe): %sA. PASSED\n',nb); else fprintf('Msg(coo_read_exe): %sA. FAILED\n',nb); end
max_difference = norm(A - B, 'inf'), if max_difference == 0, fprintf('Msg(coo_read_exe): %sB. PASSED\n',nb); else fprintf('Msg(coo_read_exe): %sB. FAILED\n',nb); end

