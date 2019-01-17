function sparse_invariants_exe()

% watch out: sparse invariants is for column vectors (v can have multiple columns)
sparse_invariants_example(3, [1; 3], [1; 2], 1, '1');
sparse_invariants_example(10, [1; 3; 4; 7; 6; 8], [1; 2; 3; 4; 5; 6], 0, '1');
sparse_invariants_example(1, [1; 3], [1; 2], 0, '1');
sparse_invariants_example(4, [1; 3], [1; 2; 3], 0, '1');
sparse_invariants_example(3, [1; 2; 3], [1, 0; 2,-1; 3,-5], 1, '1');



function sparse_invariants_example(n, i, v, okornot, nb)
check = sparse_invariants(n, i, v);
if check == okornot, fprintf('Msg(sparse_invariants_example): %s sparse_invariants PASSED\n',nb); else fprintf('Msg(sparse_invariants_example): %s sparse invariants FAILED\n',nb); end
