function ccs2mat_exe()
clear, clc

ccs2mat_example([],'1');
ccs2mat_example([1,0, 0,0,0; 3,0,-5,0,0; 0,0, 7,0,0; 0,0, 0,8,0],'2');
ccs2mat_example(toeplitz_2d_2wn(3,[-1;-2;-3;-4;-5]),'3');

function ccs2mat_example(A, nb)
full(A)
[n, m, ccs, i, v] = mat2ccs(A);
B = ccs2mat(n, m, ccs, i, v);
if norm(A - B, 'fro') > 0, fprintf('Msg(ccs2mat_exe): %s. FAILED\n',nb); else fprintf('Msg(ccs2mat_exe): %s. PASSED\n',nb); end;
