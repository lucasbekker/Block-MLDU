function ccs_scale_sparse_exe()
ccs_scale_sparse_example([], [], [], '1');
ccs_scale_sparse_example(magic(4), [], [], '2');
ccs_scale_sparse_example(magic(4), [1; 3], [2; 4], '3');

function ccs_scale_sparse_example(A, I, D, nb)
[n, m, ccs, i, v] = mat2ccs(A);
fprintf('Msg(ccs_mult_exe): Start test %s.\n',nb);
full_A = full(A)
full_I_tp = full(I)'
full_D_tp = full(D)'
n2 = m
i2 = I;
v2 = D;
d = ones(m,1); d(I) = D; d = diag(d);
full_d_scale = full(d)
[n, m, ccs, i, v] = ccs_scale_sparse(n, m, ccs, i, v, n2, i2, v2);
AD = ccs2mat(n, m, ccs, i, v); full_AD = full(AD)
if norm(AD-A*d, 'fro') > 0, fprintf('Msg(ccs_mult_exe): %s. FAILED\n',nb); else fprintf('Msg(ccs_mult_exe): %s. PASSED\n',nb); end;
