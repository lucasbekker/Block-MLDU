function coo2css_exe()
% examples
clear, clc
A = []
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 0, 0, [1], [], sparse(0,2), '1')

A = sparse(3,0)
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 3, 0, [1], [], sparse(0,2), '2')

A = zeros(0,3)
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 0, 3, [1], [], sparse(0,2), '3')

% matrix with mostly zero upper triangular part crs of that part only stores rows 1 & 2 whereas ccs stores columns 1:6 ...
A = tril(magic(7),-1); A(2,2) = 1
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 7, 7, [1;7;13;17;20;22;23], [2;3;4;5;6;7;2;3;4;5;6;7;4;5;6;7;5;6;7;6;7;7], [[38;46;5;13;21;22;1;6;14;15;23;31;16;24;32;40;33;41;49;43;2;11],[0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]], '4')

% transpose of previous case, now crs will store rows 1:7 and ccs will be empty ...
A = tril(magic(7),-1); A(2,2) = 1; A = A'
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 7, 7, [1;7;13;17;20;22;23], [2;3;4;5;6;7;2;3;4;5;6;7;4;5;6;7;5;6;7;6;7;7], [[0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0],[38;46;5;13;21;22;1;6;14;15;23;31;16;24;32;40;33;41;49;43;2;11]], '5')

% square sparse toeplitz matrix
A = toeplitz([1, 0, 3, 0, 5], [1, -1, 0, 0, 7])
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 5, 5, [1;5;8;11;13;14], [1;2;3;5;2;3;4;3;4;5;4;5;5], [1,0,3,5,1,0,3,1,0,3,1,0,1;1,-1,0,7, 1,-1,0,1,-1,0,1,-1,1]', '6')

% non square n x m with n < m matrix
A = [1, 2, 3, 4; 5, 0, 6, 7]
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 2, 4, [1;5;7], [1;2;3;4;3;4], [[1;5;0;0;0;0],[1;2;3;4;6;7]], '7')

% non square n x m with n > m matrix
A = [1, 2, 3, 4; 5, 0, 6, 7]'
[n, m, i, j, v] = mat2coo(A);
[n, m, css, i, lu] = coo2css(n, m, i, j, v); n, m, css', i', lu'
coo2css_invs(n, m, css, i, lu, 4, 2, [1;5;7], [1;2;3;4;3;4], [[1;2;3;4;6;7],[1;5;0;0;0;0]], '8')

function coo2css_invs(n, m, css, i, lu, nok, mok, cssok, iok, luok, nb)
if ccs_equal(n, m, css, i, lu, nok, mok, cssok, iok, luok) , fprintf('Msg(coo2css_exe): %s. PASSED\n',nb); else fprintf('Msg(coo2css_exe): %s. FAILED\n',nb); end
