function block_trim_exe()
block_trim_example([], 0, [], '0');
block_trim_example([], 3, [1, 1, 1], '1');
block_trim_example([1, 2, 1, 1, 1], 8, [1, 2, 1, 1, 1, 1, 1], '2');
block_trim_example([4, 2, 2], 8 , [4, 2, 2], '3');
block_trim_example([1, 3, 4, 2, 3], 7, [1, 3, 3], '4');
block_trim_example([1, 3, 4, 2, 3], 8, [1, 3, 4], '5');
block_trim_example([1, 3, 4, 2, 3], 9, [1, 3, 4, 1], '6');
block_trim_example([1, 3, 4, 2, 3], 11, [1, 3, 4, 2, 1], '7');

function s = block_trim_example(s, mn, sok, nb)
fprintf('In s: \n'); s
fprintf('In mn: \n'); mn
s = diff(block_trim(cumsum([1, s]), mn));
fprintf('Out s: \n'); s
if (length(s) == length(sok)) && (norm(s - sok, 'inf') == 0), fprintf('Msg(block_trim_exe): %s PASSED\n',nb); else fprintf('Msg(block_trim_exe): %s FAILED\n',nb); end

