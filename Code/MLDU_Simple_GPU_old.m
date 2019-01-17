function [ L, D, U ] = MLDU_Simple_GPU_old( A, s )

% Creating indices structure
j = 0;
[m,n] = size(A);
i_j = struct('i',[],'j',[]);
index = struct('M',i_j,'L',i_j,'D',i_j,'U',i_j);

% Test for sparse
if issparse(A)
    A = full(A);
end

% Convert to gpuArray
A     = gpuArray(A);

% Preallocate output matrices
L = zeros(m,n,'gpuArray');
D = zeros(m,n,'gpuArray');
U = zeros(m,n,'gpuArray');

% Loop over the blocks
for i = 1:length(s)
   
    update_index_nested(i);
     
    L(index.L.i,index.L.j) = A(index.L.i,index.L.j);
    D(index.D.i,index.D.j) = A(index.D.i,index.D.j);
    U(index.U.i,index.U.j) = A(index.U.i,index.U.j);
    
    schur_complement =  A(index.L.i,index.L.j)* ...
                       (A(index.D.i,index.D.j)\ ...
                        A(index.U.i,index.U.j));
    
    A(index.M.i,index.M.j) = A(index.M.i,index.M.j) - schur_complement;
        
end

% Nested index function
function [ ] = update_index_nested(i)
   
    index.D.i = (1:s(i)) + j;
    index.D.j = index.D.i;
    index.L.i = (index.D.i(end) + 1):m;
    index.L.j = index.D.j;
    index.U.i = index.D.i;
    index.U.j = (index.D.j(end) + 1):n;
    index.M.i = index.L.i;
    index.M.j = index.U.j;
    
    j = s(i) + j;
    
end

end








