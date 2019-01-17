function mat2coo_exe()

% examples
mat2coo_invs([], [], [], [], '1')
mat2coo_invs(sparse(3,0), [], [], [], '2')
mat2coo_invs((1:10)', [1:10]', ones(10,1), [1:10]', '3')
mat2coo_invs((1:10), ones(10,1), [1:10]', [1:10]', '3')

function mat2coo_invs(A, iok, jok, vok, nb)
[n, m, i, j, v] = mat2coo(A); n, m, i', j', v'
if abs(n-size(A,1)) + abs(m-size(A,2)) + max(abs(i - iok)) + max(abs(j-jok)) + max(max(abs(v - vok))) > 0 , fprintf('Msg(mat2coo_exe): %s. FAILED\n',nb); else fprintf('Msg(mat2coo_exe): %s. PASSED\n',nb); end
