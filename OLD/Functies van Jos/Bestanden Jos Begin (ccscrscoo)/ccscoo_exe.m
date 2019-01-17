fprintf('Basic matlab sparse matrix/vector operations and expectations (PASSED means: expectations met):\n');
fprintf('[] array differs from scalar 0:\n');
if [] == 0, fprintf('Msg(ccscoo): FAILED\n'); else fprintf('Msg(ccscoo): PASSED\n'); end
fprintf('[0] array equals to scalar 0:\n');
if [0] == 0, fprintf('Msg(ccscoo): PASSED\n'); else fprintf('Msg(ccscoo): FAILED\n'); end
fprintf('[1,2] array equals itself:\n');
if [1,2] == [1,2], fprintf('Msg(ccscoo): PASSED\n'); else fprintf('Msg(ccscoo): FAILED\n'); end
fprintf('1+[] and [] are both 0x0 matrices, but they turn out to be NOT EQUAL ... (octave3.8.1):\n');
[]   % -> [](0x0)
1+[] % -> [](0x0)
if 1+[] == [] , fprintf('Msg(ccscoo): FAILED\n'); else fprintf('Msg(ccscoo): PASSED\n'); end % this is a miracle ?
fprintf('c - A = [c - a_{ij}]_{ij}:\n');
if [1]-[1;1;1;1] == [0;0;0;0] , fprintf('Msg(ccscoo): PASSED\n'); else fprintf('Msg(ccscoo): FAILED\n'); end % this is a miracle ?

fprintf('Transformations between matlab (sparse) matrix and sparse formats:\n');
fprintf('Transformations matlab (sparse) --> sparse formats:\n');
mat2coo_exe();
mat2ccs_exe();
mat2css_exe();
fprintf('Transformations sparse formats  --> matlab (sparse):\n');
ccs2coo_exe();
coo2ccs_exe();
coo2css_exe();
coo2mat_exe();
ccs2mat_exe();
css2mat_exe();

fprintf('Sparse vector operations:\n');
sparse_zeros_exe();
sparse_equal_exe();
sparse_fxy_exe();
sparse_join_exe();
sparse_find_exe();
block_trim_exe();

fprintf('Sparse vector based Matrix operations + -:\n');
ccs_zeros_exe();
ccs_equal_exe();
ccs_get_sparse_exe();
ccs_set_sparse_exe();
ccs_fxy_exe();
ccs_add_exe();
ccs_mult_exe();
fprintf('Sparse vector based Matrix operations: x = T inv B:\n');
cxs_T_inv_cxs_B_exe();
cys_T_inv_cys_B_exe();
fprintf('Sparse vector based Matrix operations: X = l dinv u factorization:\n');
css_ldu_exe()
css_mldu_exe()

fprintf('Input/Output -- & destroy created files:\n');
mat_write_exe();
unlink('b.dat');
coo_write_exe();
unlink('A.dat');
coo_read_exe();
unlink('coo.data');
ccs_write_exe();
unlink('X.dat');
ccs_read_exe();
unlink('ccsc.data');
unlink('ccs.data');
