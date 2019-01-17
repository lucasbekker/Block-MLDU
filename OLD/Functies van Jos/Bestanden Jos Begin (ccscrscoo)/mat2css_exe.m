function mat2css_exe()
% examples
clear, clc
mat2css_invs([], [1], [], sparse(0,2), '1')
mat2css_invs(sparse(3,0), [1], [], sparse(0,2), '2')
mat2css_invs(zeros(0,3), [1], [], sparse(0,2), '3')
mat2css_invs(zeros(3,0), [1], [], sparse(0,2), '4')

fprintf('matrix with mostly zero upper triangular part crs of that part only stores rows 1 & 2 whereas ccs stores columns 1:6 ...:\n');
A = tril(magic(7),-1), A(2,2) = 1
mat2css_invs(A, [1;7;13;17;20;22;23], [2;3;4;5;6;7;2;3;4;5;6;7;4;5;6;7;5;6;7;6;7;7], [[38;46;5;13;21;22;1;6;14;15;23;31;16;24;32;40;33;41;49;43;2;11],[0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]], '5')

fprintf('transpose of previous case, now crs will store rows 1:7 and ccs will be empty ...:\n');
A = tril(magic(7),-1); A(2,2) = 1; A = A'
mat2css_invs(A, [1;7;13;17;20;22;23], [2;3;4;5;6;7;2;3;4;5;6;7;4;5;6;7;5;6;7;6;7;7], [[0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0],[38;46;5;13;21;22;1;6;14;15;23;31;16;24;32;40;33;41;49;43;2;11]], '6')

fprintf('square sparse toeplitz matrix:\n');
mat2css_invs(toeplitz([1, 0, 3, 0, 5], [1, -1, 0, 0, 7]), [1;5;8;11;13;14], [1;2;3;5;2;3;4;3;4;5;4;5;5], [1,0,3,5,1,0,3,1,0,3,1,0,1;1,-1,0,7, 1,-1,0,1,-1,0,1,-1,1]', '7')

fprintf('non square n x m with n < m matrix:\n');
mat2css_invs([1, 2, 3, 4; 5, 0, 6, 7], [1;5;7], [1;2;3;4;3;4], [[1;5;0;0;0;0],[1;2;3;4;6;7]], '8')

fprintf('non square n x m with n > m matrix:\n');
mat2css_invs([1, 2, 3, 4; 5, 0, 6, 7]', [1;5;7], [1;2;3;4;3;4], [[1;2;3;4;6;7],[1;5;0;0;0;0]], '9')

function mat2css_invs(A, cssok, iok, vok, nb)
[n, m, css, i, v] = mat2css(A); n, m, css', i', v'
if ccs_equal(n, m, css, i, v, size(A,1), size(A,2), cssok, iok, vok) , fprintf('Msg(mat2ccs_exe): %s. PASSED\n',nb); else fprintf('Msg(mat2ccs_exe): %s. FAILED\n',nb); end
