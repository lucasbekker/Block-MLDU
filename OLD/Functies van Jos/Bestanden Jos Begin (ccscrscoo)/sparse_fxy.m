% % sparse_fxy can be used to implement sparse addition, scaled addition, hadarmard multiplication, etc
%
% % TYPICAL USE: calculate z = v+w, z = v+a*w, etc.
% %              szr: store zero results
% %              szr == 0: output / store only non-zero z_i values
% %              szr == 1: output / store all  (potentially zero) z_i values
% %                         and output can contain zeros (set szr == 1)
% %              is more efficient for fxy than sparse_join followed by a fxy
% %              because only one output vector is formed in addition and if wished without zeros
%
%
function [n, i, v] = sparse_fxy(n1, i1, v1, n2, i2, v2, f, szr, echo)
if (nargin < 8) || isempty(szr), szr = 0; end;
if (nargin < 9) || isempty(echo), echo = 0; end;

sparse_invariants(n1, i1, v1);
sparse_invariants(n2, i2, v2);
if echo > 0
 fprintf('Msg(sparse_fxy): Start ... To be fxy-ed:\n');
 sparse_write(n1, i1, v1, 'v1');
 sparse_write(n2, i2, v2, 'v2');
end


[n, i, v] = sparse_zeros(max(n1, n2));

% obtain the zero element related to both arrays
zero1 = zeros(1,max(1,size(v1,2)));
zero2 = zeros(1,max(1,size(v2,2)));


ni1 = length(i1); ni2 = length(i2);
l = 1; r = 1;
while  (l <= ni1) && (r <= ni2)
 d = i1(l) - i2(r);
 if     d == 0
  e = f(v1(l,:),v2(r,:));
  if max(abs(e)) > 0 || szr > 0
   [n, i, v] = sparse_stack(n, i, v, i1(l), e);
  end
  l = l+1; r = r+1; 
 elseif d > 0
  e = f(zero1,v2(r,:));
  if max(abs(e)) > 0 || szr > 0
   [n, i, v] = sparse_stack(n, i, v, i2(r), e);
  end
  r = r+1;
 else
  e = f(v1(l,:),zero2);
  if max(abs(e)) > 0 || szr > 0
   [n, i, v] = sparse_stack(n, i, v, i1(l), e);
  end
  l = l+1;
 end
end
for k=l:ni1
 e = f(v1(k,:),zero2);
 if max(abs(e)) > 0 || szr > 0
  [n, i, v] = sparse_stack(n, i, v, i1(k), e);
 end
end
for k=r:ni2
 e = f(zero1,v2(k,:));
 if max(abs(e)) > 0 || szr > 0
  [n, i, v] = sparse_stack(n, i, v, i2(k), e);
 end
end
