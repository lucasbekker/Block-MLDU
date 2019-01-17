The executable which tests function X.m has name X_exe.m:

block_trim_exe.m      ccs_mult_exe.m          coo2css_exe.m          cxs_T_inv_cxs_B.m      sparse_find_exe.m
block_trim.m          ccs_mult.m              coo2css.m              cys_T_inv_cys_B_exe.m  sparse_find.m
ccs2coo_exe.m         ccs_read_exe.m          coo2mat_exe.m          cys_T_inv_cys_B.m      sparse_fxy_exe.m
ccs2coo.m             ccs_read.m              coo2mat.m              mat2ccs_exe.m          sparse_fxy.m
ccs2mat_exe.m         ccs_scale_sparse_exe.m  coo_read_exe.m         mat2ccs.m              sparse_invariants_exe.m
ccs2mat.m             ccs_scale_sparse.m      coo_read.m             mat2coo_exe.m          sparse_invariants.m
ccs_add_exe.m         ccs_set_sparse_exe.m    coo_write_exe.m        mat2coo.m              sparse_join_exe.m
ccs_add.m             ccs_set_sparse.m        coo_write.m            mat2css_exe.m          sparse_join.m
ccscoo_exe.m          ccs_stack_sparse.m      css2mat_exe.m          mat2css.m              sparse_stack.m
ccs_copy.m            ccs_truncate.m          css2mat.m              mat_write_exe.m        sparse_write.m
ccs_equal_exe.m       ccs_write_exe.m         css_ldu_exe.m          mat_write.m            sparse_zeros_exe.m
ccs_equal.m           ccs_write.m             css_ldu_exe.ok         number_equal.m         sparse_zeros.m
ccs_fxy_exe.m         ccs_zeros_exe.m         css_ldu.m              README.txt             stokes12.coo
ccs_fxy.m             ccs_zeros.m             css_mldu_exe.m         sequence_equal.m       stokes54.coo
ccs_get_sparse_exe.m  coo2ccs_exe.m           css_mldu.m             sparse_equal_exe.m     TODO.txt
ccs_get_sparse.m      coo2ccs.m               cxs_T_inv_cxs_B_exe.m  sparse_equal.m

THE ROUTINE to transform into C/C++ is: css_mldu.m
All other routines on ccs matrices such as:
 ...
 ccs_add.m   (A+B)
 ccs_copy.m  (A = B)
 ccs_equal.m
 ccs_mult.m
 ...
are test the datastructures











implementation issues
---------------------
we use matlab vectors i = [] for indices (which default to 0 x 0 matrices)
however, we should have used empy matlab VECTORS i = sparse(0,1) for indices ...
this creates trouble in array operations such as (abs(v) > 0) & (i > j)
for the empty matrix ... We patched in in css2mat by adding a switch: if length(i)> 0 else end...
in coo2css there is the reset to i = [] after i = [] and adding to empty slots yields i = 2 x 0 ...
(instead of the desired 0 x 1) -- but this is a matlab convention which just turns out to be unhandy
for our specific interest ...




matlab/octave bugs
------------------
* octave3.8.1> v = []; norm(v) -> ans =   7.4817e-317
   (norm(v,'inf') -> 0, max(abs(v)) -> 0 are ok)
* matlab claims that "L = length(X) returns the length of the largest array dimension in X", which is not true.
   Since matlab in addition claims: "The length of an empty array is zero." (empty array is 0 x m x ... or n x 0 x ... etc) 
    As a result, the length of a column vector can DIFFER from size(X,1) (it is zero if size(X,2) = 0 ...)
     and one MUST use length(X) instead of size(X,1) to get the proper length of a vector


matlab commands called and matlab implementation properties relied on
---------------------------------------------------------------------
- sparse(n, m): 
  * empty sparse matrix + c -> empty sparse matrix
  * empty sparse matrix * c -> empty sparse matrix

- index-subset selection i(begin:1:end):
  * v(k>l) -> empty sparse vector, also if v is empty sparse matrix (vector)


matlab sparse vector/matrix toolbox
-----------------------------------
implements sparse vector/matrix operations in matlab
is excruciatingly slow due to amount of indexing

reason to exist:
- for educational purposes
- for prototypes -- before writing c/c++/fortran code

sparse vector definition
------------------------
is a triple(i, v, n) such that invariants:
0. n >= 0
1. 0 < i_k < i_{k+1} for k=1:|i|-1 [lower bound not checked]
2. |v| == |i|
3. i_k > n ==> v_k = 0

?. i_k == i_{k+1} is (see 1.) forbidden. Is this really necessary?
!. i_k > n ==> v_k = 0 is needed for css storage when n<>m

Remark: A sparse vector is (also) used as a sparse diagonal matrix
        where missing entries have value 1 (???????????????)
        in ccs_scale_sparse()

sparse vector operations
------------------------
sparse_zeros: create a zero sparse vector of given length n
sparse_equal
sparse_invariants: checks all above sparse vector invariants
sparse_stack
sparse_join: joins two sparse column vectors/matrices into one sparse column matrix
sparse_fxy:
sparse_write:
sparse_map: to be implemented:
  is trivial if f(0) is 0, then (i, v, n) -> (i, f(v), n).
  however, if f(0) <> 0, then sparse vector should become full ...
sparse vectors are used to build sparsely represented matrices:

sparse matrix formats internal
------------------------------
my coo
------
[n, m, i, j, v]
i,j,v (column) vectors or empty vectors: SHOULD HAVE BEEN ROW VECTORS FROM THE OLD xroutine() MAPPING POINT OF VIEW!!!!!!!!!!
with invariants:
? n, m integer >= 0 (?) not checked
- |i| = |j| = |v|
- i_k > n ==> v_k = 0 (weaker than max demand i_k <= n, matlab sparse() automatically deletes the zero entries in conversion to matlab-sparse matrix)
- j_k > n ==> v_k = 0 (weaker than max demand i_k <= n, matlab sparse() automatically deletes the zero entries in conversion to matlab-sparse matrix)
? (i_k, j_k) <> (i_l, j_l) for all k <> l (different in Euclidean distant (DO WE USE THIS? -- matlab sparse() can handle double index pairs, simply adds up the related values))
? ((j_k == j_{k-1} & i_k > i_{k-1}) || (j_k > j_{k-1})) for all k=1:|i|-1 (DO WE USE THIS? -- THIS INV is what matlab output [i,j,v] =find(A) satisfies)

my cxs (ccs/crs, x =[c|r])
----------
compressed column storage:
[n, m, cxs, i, v]
stores n x m matrix A in ccs format exor
stores m x n matrix A^T in css format
[<==> stores n x m matrix A in crs format]
with invariants:
- |i| = |v|
- |i| = ccs(end) - 1 (<=> ccs(|ccs|) = |i| + 1)
- |ccs|-1 <= m;
- [n, i(ccs(k)):i(ccs(k+1)-1), v(ccs(k)):v(ccs(k+1)-1)] is a sparse vector for all k = 1:|ccs|-1
- [n, m, ccs, i, lu(:,1)] satisfies css invariants (and forms a lower triangular matrix)
- [m, n, ccs, i, lu(:,1)] satisfies css invariants (and forms a upper triangular matrix)

 EXAMPLE: n = 3, m = 4, |ccs|-1 = 2 is the ``matrix'' with potential nonzeros at positions (*)

                                 ****
                                 ****
                                 **

 and is helpfull for matrices related to constrained problems, lower right block needs not to be stored


 REMARK: no forced entry storage at certain position such as the diagonal (i,i) location
 REMARK: Let v be Natural^{m x 1} vector. Define w := cumsum([1;v]). Then v = diff(w) -- matlab diff and cumsum command.



my css
------
purpose of existing: index-searching free Gaussian Elimination implementation possible if matrix
is stored in css format.
compressed symmetric storage (css):
[n, m, css, i, lu]
with invariants:
1. i_k <= max(n, m) for all k=1:css(end)-1     -- needs to be able to store all non-zeros of an n <> m matrix
3. i(css(k)), ..., i(ccs(k+1)-1) >= k for all k=1,...,|css|-1 (only store lower triangular part)
4. [n, m, css, i, lu(:,1)] L+D    lower triangular (with diagonal D and strictly lower triangular part L)
                          in css format satisfies ccs invariants
5. [m, n, css, i, lu(:,2)] (U+D)' upper triangular (with diagonal D and strictly upper triangular part U)
                          in ccs format (<==> U upper triangular in crs format) satisfies ccs invariants
2. c := |css|-1 <= min(n, m) -- follows from invariants 4. and 5.

- if n <> m (ONLY) zero entries can be stored outside n x m block
- stored scalar-diagonal D entries are doubly stored (right now) i.e., l_kk = u_kk = d_kk if entry kk is stored
- by construction, potential diagonal k x k blocks (k >= 1) do contain entries in both L+D and U+D


 REMARK: Let n = 3, m = 4, c = 2. Stored is a matrix with potential nonzeros at positions (*)

                                 ****
                                 ****
                                 **

         and is helpfull for matrices related to constrained problems.

 REMARK: css stores extra zeros, at least s extra zeros (for diagonal entries a_11 ... a_ss)
         need not store these zeros, but simplifies current implementation
 REMARK: no forced entry storage at certain position such as the diagonal (i,i) location
 REMARK: i(css(k)) = k <==> diagonal entry a_kk is stored at both lu(i(css(k)), 1) and lu(i(css(k)), 2)
          (if no entry is stored for row/col k then i(css(k)) == i(css(l)) for some row/col l > k
           and since l > k one would find (if l is stored) i(css(k)) = i(css(l)) >= l > k)


disk io convention
------------------
write column AND row vectors in V^n with V in {integer, real, complex} as row vectors to disk. This way only one 'newline' or 'newline\return'
combination is needed per vector (if the vector is written as a column vector to disk, each entry has to be preceeded by a 'newline'.

sparse matrix conversion and to/from disk routines
--------------------------------------------------
X2Y: matlab function which transforms format X into format Y [formatY] = X2Y(formatX)
X_write: matlab function which writes format X into format Y on disk []=X2Y(formatX, 'name) -> formatY on disk in file 'name'
X_read: matlab function which reads  format X from disk from file 'name' into format [formatY]=X2Y('name')


sparse matrix formats to and from disk
--------------------------------------
mat: R + W
coo: R + W
ccs: R + W


Sparse matrix formats: definition
---------------------------------

ccs format sparse vector
------------------------
[n, m, ccs, i, v]:
n x m is size matrix
ccs(j), j=1,...,(end-1)=m: amount of entries in column j
i(ccs(j):ccs(j+1)-1): row entries i in column j
v(ccs(j):ccs(j+1)-1): values v_ij related entries ij in column j
such that triples (i(ccs(j):ccs(j+1)-1), v(ccs(j):ccs(j+1)-1), n) satisfy sparse vector definition for all j=1,...,m

my crs
------
compressed row       storage (crs): is identical to ccs storage for A'
so far: [n, m, crs, i, v]
FUTURE? [m, n, crs, i, v] IS BETTER BECAUSE set/get on row/column can then give back sparse sequence of length first arg


sparse matrix operations: ccs operations
----------------------------------------
ccs_zeros
ccs_equal
ccs_copy: use to copy a matrix for in-situ factorizations
ccs_stack: stack (append) sparse column to ccs matrix
ccs_set_sparse
ccs_get_sparse
ccs_fxy: implements ccs_add
ccs_mult
ccs_scale_sparse: not implemented is trivial (v := v*c)


Sparse LU implementations
-------------------------

UMFPACK uses multifrontal LU approach (calculates entire Schur complement updates)
SuperLU: super nodal approach (serial mode)
only ma77 seems to take inputted pivot pairs and stick to them

HSL drivers
-----------

hsl drivers run ma*7 with * in [5, 7, 9].

All read X.dat (matrix) and b.dat (rhs vector(s))

and return x.dat (solution).

Ma97 driver returns also P.dat (permutation)

57:
---

77:
---
X.dat in ccsc (compressed column storage chunck, see ".ds" example ./examples/Fortran) format
       (see HSL77 2013, Sec4/MA77_input_vars)

97:
---
X.dat: must contain 'tril(X)' in ccs (compressed column storage, see ".ds" example ./examples/Fortran) format
       (see HSL97 2013, Sec2.6.3)
REMARK: mau driver reads full ccs format 'n,m,nne', whereas native hsl example ".ds"
         only reads castrated ccs format variant 'n,nne'
          this castrated format is possible because 97 assumes X = X^T ==> n == m
REMARK: take care that 'X.dat' only contains ccs of 'tril(X)' -- if only to avoid the 'out of range' warning


