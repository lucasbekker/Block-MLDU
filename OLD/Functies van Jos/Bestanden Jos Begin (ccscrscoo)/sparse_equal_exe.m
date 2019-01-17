function sparse_equal_exe()
fprintf('Msg(sparse_equal): Sparse vectors v1 and v2 are equal:\n');
sparse_equal_invs(4, [1; 3], [-1; 8], 4, [1; 3], [-1; 8], '1')
fprintf('Msg(sparse_equal): Sparse vectors v1 and v2 are different:\n');
sparse_equal_invs(4, [1; 8], [-1; 8], 4, [1; 3], [-1; 8], '2',1)

function sparse_equal_invs(n1, i1, v1, n2, i2, v2, nb, fail)
if nargin < 8 || isempty(fail), fail = 0; end;
if fail == 0
 if  sparse_equal(n1, i1, v1, n2, i2, v2) , fprintf('Msg(sparse_equal_exe): %s. PASSED\n',nb); else fprintf('Msg(sparse_equal_exe): %s. FAILED\n',nb); end
else
 if ~sparse_equal(n1, i1, v1, n2, i2, v2) , fprintf('Msg(sparse_equal_exe): %s. PASSED\n',nb); else fprintf('Msg(sparse_equal_exe): %s. FAILED\n',nb); end
end
