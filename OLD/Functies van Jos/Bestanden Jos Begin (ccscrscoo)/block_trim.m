% input + output IN cumsum format!
function s = block_trim(s, mn)
% for square matrices, s(end) = n + 1
if     s(end) < mn + 1
 % too short, add 1x1 blocks
 s = cumsum([1, diff(s), ones(1, mn + 1 - s(end))]);
elseif s(end) > mn + 1
 % too long, truncate (inside) the first block which exceeds 
 f = sparse_find(mn + 1, s)
 s = cumsum([1, diff(s(1:f)), mn + 1 - s(f)]);
end
