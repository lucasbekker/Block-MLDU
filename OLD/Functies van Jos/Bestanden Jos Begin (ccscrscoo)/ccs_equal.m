% Tests only whether storage fields are exactly the same
% Test can be used for ccs AND crs (n and m do not play different roles)
%
% Tests NOT whether represented matrices might be the same (could be 
% since ccs can store (or not) empty columns).
function eql = ccs_equal(n1, m1, ccs1, i1, v1, n2, m2, ccs2, i2, v2, epsilon, echo)
if (nargin < 11) || isempty(epsilon), epsilon = 2^-44; end;
if (nargin < 12) || isempty(echo), echo = 0; end;
% last two checks should be part of ccs_invariants!
eql = (number_equal(n1, n2, 0, echo) && number_equal(m1, m2, 0, echo+10) && sequence_equal(ccs1, ccs2, 0, echo) && sequence_equal(i1, i2, 0, echo+100) && sequence_equal(v1, v2, [], echo+200) && number_equal(length(i1)*size(i1, 1), length(i1)*size(v1, 1), 0, echo+1000) && number_equal(length(i1)*size(i1, 1), length(i2)*size(v2, 1), 0, echo+1001));
