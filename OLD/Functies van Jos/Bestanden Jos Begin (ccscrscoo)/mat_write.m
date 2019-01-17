% write n x m matrix full format (all entries) into file
% matlab command dlmwrite always writes the full matrix,
% even if the matrix is sparse. matlab command disp
% would print a sparse matrix in [i, j, v] format
%
function mat_write(A, fname)
if (nargin < 2) || (isempty(fname))
 disp(full(A));
else
 dlmwrite(fname,A,'delimiter',' ','precision','%22.16e');
end
