function [n, m, i, j, v] = mat2coo(A)
[n, m] = size(A);
[i,j,v] = find(A);
% matlab returns vectors i, j, v of size 0x1 if there are all entries are zero -- this could be called correct (indicates a column vector of length zero)
% testing against the empty vector we do with testing against [] (but that turns out to be a 0x0 vector) would go wrong -- fortran code should define []-vector as 0x1 vector
% and all is then fine again.
if     length(i) == 0
 i = []; j = []; v = [];
% matlab return column vectors i, j, v except for a matrix which consists of one single row ...
% our coo format relies on column vectors:
elseif size(i,2) > size(i,1)
 i = i'; j = j'; v = v';
end;
% which are matlab 0x0 empty vectors (and not for instance 0x1 empty vectors ...)
if length(A) == 0, i = []; j = []; v = []; end
