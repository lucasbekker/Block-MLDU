function sparse_find_exe()
fprintf('sparse empty vectors:\n');
sparse_find_invs(1, [], 0, '1');

fprintf('sparse 3 x 1 vectors:\n');
sparse_find_invs(15, [3;8;11], 3, '2');

fprintf('sparse 4 x 1 vectors:\n');
sparse_find_invs(15, [3;8;11;14], 4, '3');

fprintf('sparse 5 x 1 vectors, double entries i(p) == k = %d:\n',15);
sparse_find_invs(15, [3;8;11;14;14], 5, '4');

fprintf('sparse 5 x 1 vectors, all smaller than k = %d:\n',16);
sparse_find_invs(16, [3;8;11;14;14], 5, '5');

fprintf('sparse 5 x 1 vectors:\n');
sparse_find_invs(14, [3;8;11;14;14], 3, '6');
sparse_find_invs(14, [14;14], 0, '7');

fprintf('sparse 5 x 1 vectors, all larger than k = %d:\n',3);
sparse_find_invs(3, [3;8;11;14;14], 0, '7')

fprintf('sparse 5 x 1 vectors:\n');
sparse_find_invs(4, [3;8;11;14;14], 1, '8');


function sparse_find_invs(k, i, lok, nb)
l = sparse_find(k, i)
if abs(l-lok) == 0 , fprintf('Msg(sparse_find_exe): %s. PASSED\n',nb); else fprintf('Msg(sparse_find_exe): %s. FAILED\n',nb); end

