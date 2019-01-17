function ccs_set_sparse_exe()
[n, m, ccs, i, v] = ccs_zeros(4, 3);
[n, m, ccs, i, v] = ccs_set_sparse_example(n, m, ccs, i, v, 2, 4, [2;4], [-1;-5],4, 3, [1;1;3],[2;4],[-1;-5], '1');
[n, m, ccs, i, v] = ccs_set_sparse_example(n, m, ccs, i, v, 1, 5, [1;3;5], [-7;-2;-4], 5, 3, [1;4;6], [ 1;3;5;2;4], [-7;-2;-4;-1;-5], '2');
[n, m, ccs, i, v] = ccs_set_sparse_example(n, m, ccs, i, v, 4, 5, [1;2;5], [-1;-2;-5], 5, 4, [1;4;6;6;9], [1;3;5;2;4;1;2;5], [-7;-2;-4;-1;-5;-1;-2;-5], '3');
[n, m, ccs, i, v] = ccs_set_sparse_example(n, m, ccs, i, v, 1, 4, [], [], 5, 4, [ 1;1;3;3;6], [ 2;4;1;2;5], [-1;-5;-1;-2;-5],'4');

function [n, m, ccs, i, v] = ccs_set_sparse_example(n, m, ccs, i, v, k, nk, ik, vk, nok, mok, ccsok, iok, vok, nb)
n, m, ccs_tp = ccs', i_tp = i', v_tp = v', k, nk, ik_tp = ik', vk_tp = vk'
[n, m, ccs, i, v] = ccs_set_sparse(n, m, ccs, i, v, k, nk, ik, vk);
n, m, ccs_tp = ccs', i_tp = i', v_tp = v', nok, mok, ccsok_tp = ccsok', iok_tp = iok', vok_tp = vok'
if ccs_equal(n, m, ccs, i, v, nok, mok, ccsok, iok, vok) > 0 , fprintf('Msg(ccs_set_sparse_exe): %s. PASSED\n',nb); else fprintf('Msg(ccs_set_sparse_exe): %s. FAILED\n',nb); end
