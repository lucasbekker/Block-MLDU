function eql = number_equal(n1, n2, epsilon, echo)
if (nargin < 3) || isempty(epsilon), epsilon = 2^-44; end;
if (nargin < 4) || isempty(echo), echo = 0; end;
eql = (norm(n1 - n2, 'inf') <= epsilon);
if (echo > 0 || ~eql)
 fprintf('Msg(number_equal-%d): integers are equal: %d (0: NO; 1: YES)\n', echo, eql);
 fprintf('Msg(number_equal-%d): integer %d = n1 <> n2 = %d\n', echo, n1, n2);
end
