function ccs_equal_exe()

A = [0, 0, 0, 0; 0, 2, 3, 0; 0, 0, 0, 0; 0, 0, 0, 0]
fprintf('The above matrix in ccs format: \n');
n1 = 4
m1 = 4
ccs1 =[1;1;2;3]
i1 = [2;2]
v1 = [2;3]
n2 = 4;
m2 = 4;
ccs2 =[1;1;2;3];
i2 = [2;2];
v2 = [2;3];
fprintf('Msg(ccs_equal): Matrix A and B are ccs equal:\n');
ccs_equal_example(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, '1')

B = [0, 7, 0, 0; 0, 2, 3, 0; 0, 0, 0, 0; 0, 0, 0, 0]
fprintf('The above matrix in ccs format: \n');
n2 = 4
m2 = 4
ccs2 = [1;1;3;4]
i2 = [1;2;2]
v2 = [7;2;3]
fprintf('Msg(ccs_equal): Matrix A and B are ccs different:\n');
ccs_equal_example(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, '2',1)

function ccs_equal_example(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, nb, fail)
if nargin < 12 || isempty(fail), fail = 0; end;
if fail == 0
 if  ccs_equal(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2) , fprintf('Msg(ccs_equal_exe): %s. PASSED\n',nb); else fprintf('Msg(ccs_equal_exe): %s. FAILED\n',nb); end
else
 if ~ccs_equal(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2) , fprintf('Msg(ccs_equal_exe): %s. PASSED\n',nb); else fprintf('Msg(ccs_equal_exe): %s. FAILED\n',nb); end
end
