function mat2ccs_exe()

% an empty matrix must have ccs = [1] (by construction: ccs(end) = nnz + 1 must exist)
clear, clc
mat2ccs_invs([], [1], [], [], '1')
mat2ccs_invs(sparse(0,0), [1], [], [], '2')
mat2ccs_invs(zeros(3,4), [1], [], [], '3')
mat2ccs_invs([1, 0, -2, 0, 0, 3, 0, 0, 0, -4, 0, 0, 0, 0], [1;2;2;3;3;3;4;4;4;4;5],[1;1;1;1], [1;-2;3;-4], '4')
mat2ccs_invs([1, 0, -2, 0, 0, 3, 0, 0, 0, -4, 0, 0, 0, 0]', [1;5], [1;3;6;10], [1;-2;3;-4], '5')
% a sparse n x m matrix with an entry at position (n, m)
mat2ccs_invs([1,0, 0;...
              3,0,-5;...
              0,0, 7], [1;3;3;5], [1;2;2;3], [1;3;-5;7], '6')

% a sparse n x m matrix with a ZERO entry at position (n, m)
mat2ccs_invs([1,0, 0,0,0; 3,0,-5,0,0; 0,0, 7,0,0; 0,0, 0,8,0], [1;3;3;5;6], [1;2;2;3;4], [1;3;-5;7;8], '7')
mat2ccs_invs([1,0, 0,0,0;...
             3,0,-5,0,0;...
             0,0, 7,0,0;...
             0,0, 0,8,0], [1;3;3;5;6], [1;2;2;3;4], [1;3;-5;7;8], '8')
mat2ccs_invs(toeplitz_2d_2wn(3,[-1;-2;-3;-4;-5]), ...
             [1;4;8;11;15;20;24;27;31;34], ...
             [1;2;4;1;2;3;5;2;3;6;1;4;5;7;2;4;5;6;8;3;5;6;9;4;7;8;5;7;8;9;6;8;9], ...
             [-3;-2;-1;-4;-3;-2;-1;-4;-3;-1;-5;-3;-2;-1;-5;-4;-3;-2;-1;-5;-4;-3;-1;-5;-3;-2;-5;-4;-3;-2;-5;-4;-3], '9');

% example: example 1. from ma77
mat2ccs_invs([2,3,0,0,0;3,1,4,0,6;0,4,1,5,0;0,0,5,3,0;0,6,0,0,1], [1;3;7;10;12;14],[1;2;1;2;3;5;2;3;4;3;4;2;5],[2;3;3;1;4;6;4;1;5;5;3;6;1],'10');
% example: example 1. from ma97
mat2ccs_invs([2,1,0,0,0;1,4,1,0,1;0,1,3,2,0;0,0,2,0,0;0,1,0,0,2], [1;3;7;10;11;13],[1;2;1;2;3;5;2;3;4;3;2;5],[2;1;1;4;1;1;1;3;2;2;1;2],'11');
% example: example 2. from ma97
mat2ccs_invs([1,-3,0,1,0;-3,-5,6,0,4;0,6,0,2,0;1,0,2,3,0;0,4,0,0,1], [1;4;8;10;13;15],[1;2;4;1;2;3;5;2;4;1;3;4;2;5],[1;-3;1;-3;-5;6;4;6;2;1;2;3;4;1],'12');

function mat2ccs_invs(A, ccsok, iok, vok, nb)
[n, m, ccs, i, v] = mat2ccs(A); n, m, ccs_tp = ccs', i_tp = i', v_tp = v'
if ccs_equal(n, m, ccs, i, v, size(A,1), size(A,2), ccsok, iok, vok) , fprintf('Msg(mat2ccs_exe): %s. PASSED\n',nb); else fprintf('Msg(mat2ccs_exe): %s. FAILED\n',nb); end
