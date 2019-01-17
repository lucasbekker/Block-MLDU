function eql = sparse_equal(n1, i1, v1, n2, i2, v2, epsilon, echo)
if (nargin < 7) || isempty(epsilon) epsilon = 2^-44; end;
if (nargin < 8) || isempty(echo) echo = 0; end;
eql = number_equal(n1, n2, 0, echo) && sequence_equal(i1, i2, 0, echo) && sequence_equal(v1, v2, epsilon, echo+100);
% notification of "not equal" takes place via the "_equal" routine calls
