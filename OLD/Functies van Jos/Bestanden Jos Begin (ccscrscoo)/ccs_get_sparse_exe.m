function ccs_get_sparse_exe()
n = 8; m = 5;
% ccs which only stores m - 1 (!) columns, and one empty column (the last but one)
ccs = [3;2;3;0];
i = [1;2;4;2;3;2;3;4];
v = [1;2;3;4;5;6;7;8];
ccs = cumsum([1; ccs]);

%ccs_get_sparse_example(n, m, ccs, i, v, 0, 4, [1;2;4], [1;2;3], '0'); % SHOULD AND DOES FAIL
ccs_get_sparse_example(n, m, ccs, i, v, 1, 8, [1;2;4], [1;2;3], '1');
ccs_get_sparse_example(n, m, ccs, i, v, 2, 8, [2;3], [4;5], '2');
ccs_get_sparse_example(n, m, ccs, i, v, 3, 8, [2;3;4], [6;7;8], '3');
ccs_get_sparse_example(n, m, ccs, i, v, 4, 8, sparse(0,1), sparse(0,1), '4');
ccs_get_sparse_example(n, m, ccs, i, v, 5, 8, sparse(0,1), sparse(0,1), '5');

% ccs which only stores m - 2 (!) columns, and one empty column (the last but one)
% value field with > 1 columns
n = 3; m = 4;
ccs = [3;2];
i = [1;2;3;2;3];
v = [[1,2,3,4,5];[6, 7, 8, 9, 0]]';
ccs = cumsum([1; ccs]);
ccs_get_sparse_example(n, m, ccs, i, v, 2, 3, [2;3], [4, 9; 5, 0], '6');


function ccs_get_sparse_example(n, m, ccs, i, v, k, nok, iok, vok, nb)
n, m, ccs_tp = ccs', i_tp = i', v_tp = v'
[nk, ik, vk] = ccs_get_sparse(n, m, ccs, i, v, k);
nk, ik_tp = ik', vk_tp = vk'
if sparse_equal(nk, ik, vk, nok, iok, vok) > 0 , fprintf('Msg(ccs_get_sparse_exe): %s. PASSED\n',nb); else fprintf('Msg(ccs_get_sparse_exe): %s. FAILED\n',nb); end
