% CRS and CCS do not help with the most important issue:
% assume CCS: and we do LU with row i: which entries j are to be used? we would have to search ALL columns i+1:m
%  NOTE THAT CCS ALSO STORES 0 ENTRIES for columns with NO entries, i.e, REALLY ALL COLUMNS have to be checked with default CCS, worse, a non-singular matrix has a non-zero entry in each column, so none of the columns are empty ...
%  SO: CCS fine for L-updates (column vector updates) to create schurcomplement
%      CCS is a problem for determining the columns needed of U to do the L-update
%      (ONLY WAY AROUND IS to maintain BOTH CCS AND CRS)
%  DITTO BY SYMM of the microblock LU factorization FOR CRS.
% WAY around:
%  1. shift to preprocessor (AMD/RCM) -- does not solve ANYTHING except if all to be created zeroes are created in advance
%      could lead to many zeros stored (in case of incomplete factorization)
%  2. use both CCS and CRS and do 2 x LU calculations (L*D\inv and D\invU U-updates respectively L-updates)
%  2. Assume sort/extend more expensive then 2 x LU calculations and use Maubach's SCCS scheme [i, j, aij, aji]
%      (at the expense of storing potential zeros)
%      L-updates: only performed on lower(S,-1)
%      U-updates: only performed on upper(S,0)
%     Index extendion also be done symbolicaly, i.e., in advance
%
% Matlab quirks:
%  [i,j,v] = find(A) return column vectors i, j, v EXCEPT when A = 1 x m matrix, then i, j, v are row vectors ...
%
%
% Remarks:
%  - m == length(ccs) - 1 NEED NOT BE STORED as part of ccs data structure
%  - outputs i, v are n x ? matrices (n x 1 or n x 0) with n >= 0
%    in fortran/c implementation i, v will be n-vectors
%    so already here we avoid using the 2-nd dimension
%
%
%
function [n, m, ccs, i, v] = mat2ccs(A)
[n, m, i, j, v] = mat2coo(A);
[n, m, ccs, i, v] = coo2ccs(n, m, i, j, v);
