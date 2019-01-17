function css2mat_exe()
fprintf('Msg(css2mat_exe): Tests use mat2css!:\n');
% examples
clear, clc
% sparse banded square matrix
css2mat_example(toeplitz([1, 0, 3, 0, 5], [1, -1, 0, 0, 7]), '1');

% non square 2 x 4 with n < m matrix
css2mat_example([1, 2, 3, 4; 5, 0, 6, 7], '2');

% non square 4 x 2 with n > m matrix
css2mat_example([1, 2, 3, 4; 5, 0, 6, 7]', '3');

function css2mat_example(A, nb)
full(A)
[n, m, css, i, lu] = mat2css(A); n, m, css', i', lu'
B = css2mat(n, m, css, i, lu);
if norm(A-B,'inf') > 0, fprintf('Msg(css2mat_exe): %sA. FAILED\n',nb); else fprintf('Msg(css2mat_exe): %sA. PASSED\n',nb); end
[LA, UA] = css2mat(n, m, css, i, lu);
if norm(tril(A,-1)-LA,'inf') + norm(triu(A)-UA,'inf') > 0, fprintf('Msg(css2mat_exe): %sB. FAILED\n',nb); else fprintf('Msg(css2mat_exe): %sB. PASSED\n',nb); end
