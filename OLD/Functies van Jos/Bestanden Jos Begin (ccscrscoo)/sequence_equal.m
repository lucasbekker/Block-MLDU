function eql = sequence_equal(s1, s2, epsilon, echo)
if (nargin < 3) || isempty(epsilon) epsilon = 2^-44; end;
if (nargin < 4) || isempty(echo) echo = 0; end;
eql = (norm(size(s1) - size(s2), 'inf') == 0) && (norm(s1 - s2, 'inf') <= epsilon);
if echo > 0 || ~eql
 fprintf('Msg(sequence_equal-%d): sequences are equal: %d (0: NO; 1: YES)\n', echo, eql);
 fprintf('Msg(sequence_equal-%d): sequences inf-distance tolerance: %g):\n', echo, epsilon);
 fprintf('Msg(sequence_equal-%d): size s1 = \n', echo); disp(size(s1))
 fprintf('Msg(sequence_equal-%d): size s2 = \n', echo); disp(size(s2))
 fprintf('Msg(sequence_equal-%d): s1 = \n', echo); disp(s1)
 fprintf('Msg(sequence_equal-%d): s2 = \n', echo); disp(s2)
end
