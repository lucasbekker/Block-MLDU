function sparse_zeros_exe()
sparse_zeros_invs(0, '1');
sparse_zeros_invs(4, '2');

function sparse_zeros_invs(n, nb)
[nout, i, v] = sparse_zeros(n);
if abs(nout-n) + norm(i - [],'inf') + norm(v - [],'inf') > 0 , fprintf('Msg(sparse_zeros_exe): %s. FAILED\n',nb); else fprintf('Msg(sparse_zeros_exe): %s. PASSED\n',nb); end
