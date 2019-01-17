function sparse_fxy_exe()

% add two sparse (empty) vectors of length 0
sparse_fxy_example(0, [], [], 0, [], [], @(x,y) x + y, 0, 0, [], [], '1');
sparse_fxy_example(0, [], [], 0, [], [], @(x,y) x + y, 1, 0, [], [], '2');

% add two sparse vectors empty but of positive length
sparse_fxy_example(10, [], [], 30, [], [], @(x,y) x + y, 0, 30,[], [], '3');
sparse_fxy_example(10, [], [], 30, [], [], @(x,y) x + y, 1, 30,[], [], '4');

% add two sparse vectors only with entry at end
sparse_fxy_example(21, [21], [-7], 2, [1], [8], @(x,y) x + y, 0, 21, [1;21], [8;-7], '5');
sparse_fxy_example(21, [21], [-7], 2, [1], [8], @(x,y) x + y, 1, 21, [1;21], [8;-7], '6');

sparse_fxy_example(0, [], [], 2, 1, 1, @(x,y) x + y, 0, 2, [1], [1], '7');
sparse_fxy_example(0, [], [], 2, 1, 1, @(x,y) x + y, 1, 2, [1], [1], '8');

sparse_fxy_example(0, [], [], 10, [1;4;5;7;8], [-5;2;-2;8;7], @(x,y) x + y, 0, 10, [1;4;5;7;8], [-5;2;-2;8;7], '9');
sparse_fxy_example(0, [], [], 10, [1;4;5;7;8], [-5;2;-2;8;7], @(x,y) x + y, 1, 10, [1;4;5;7;8], [-5;2;-2;8;7], '10');
 
% create to sparse vectors V and W:
A = toeplitz_2d_2wn(3,[-1;-2;-3;-4;-5]); full(A)
[i4, j4, v4] = find(A(:,4))
[i7, j7, v7] = find(A(:,7))
n = size(A,1);
sparse_fxy_example(n, i4, v4, n, i7, v7, @(x,y) x - y, 0, 9, [1;4;5;7;8], [-5;2;-2;2;2], '11');
sparse_fxy_example(n, i4, v4, n, i7, v7, @(x,y) x - y, 1, 9, [1;4;5;7;8], [-5;2;-2;2;2], '12');
 
sparse_fxy_example(15, [1;3;5],[-1;-7;12], 20, [2; 3; 4], [1; 2; -4], @(x,y) x + 2*y, 0, 20, [1;2;3;4;5], [-1;2;-3;-8;12], '13');
sparse_fxy_example(15, [1;3;5],[-1;-7;12], 20, [2; 3; 4], [1; 2; -4], @(x,y) x + 2*y, 1, 20, [1;2;3;4;5], [-1;2;-3;-8;12], '14');

% both removal of zero elements and non-removal should pass this test
sparse_fxy_example(1, [1], [-1], 1, [1], [1], @(x,y) x + y, 0, 1, [], [], '15');
sparse_fxy_example(1, [1], [-1], 1, [1], [1], @(x,y) x + y, 1, 1, [1], [0], '16');
sparse_fxy_example(10, [1; 3; 4; 6; 7], [0; -3; -4; -5; -6], 8, [3; 4; 5; 6; 7;8], [3;4; 0;5; 6;0], @(x,y) x + y, 1, 10, [1; 3; 4; 5; 6; 7; 8], [0; 0; 0; 0; 0; 0; 0], '17');
sparse_fxy_example(10, [1; 3; 4; 6; 7], [0; -3; -4; -5; -6], 8, [3; 4; 5; 6; 7;8], [3;4; 0;5; 6;0], @(x,y) x + y, 0, 10, [], [], '18');
sparse_fxy_example(10, [1; 3; 4; 6; 7], [-2; -3; -4; -5; -6], 8, [1; 2; 3; 4; 7; 8], [-11,-21; -12,-22; -13,-23; -14,-24;-15,-25; -16,-26], @(x,y) x + y, 1, 10, [1; 2; 3; 4; 6; 7; 8], [-13,-23;-12,-22;-16,-26;-18,-28;-5,-5;-21,-31;-16,-26], '19');

function sparse_fxy_example(n1, i1, v1, n2, i2, v2, f, szr, nok, iok, vok, nb)
[n, i, v] = sparse_fxy(n1, i1, v1, n2, i2, v2, f, szr); n, i_tp = i', v_tp = v', nok, iok_tp = iok', vok_tp = vok'
%if abs(n-nok) + norm(i - iok,'inf') + norm(v - vok,'inf') > 0 , fprintf('Msg(sparse_fxy_exe): %s. FAILED\n',nb); else fprintf('Msg(sparse_fxy_exe): %s. PASSED\n',nb); end
if sparse_equal(n, i, v, nok, iok, vok) , fprintf('Msg(sparse_fxy_exe): %s. PASSED\n',nb); else fprintf('Msg(sparse_fxy_exe): %s. FAILED\n',nb); end
