function [ L, D, U ] = MLDU_Simple_GPU( A, s )

% Creating indices structure
j = 0;
[m,n] = size(A);
s_l = length(s);
i_j = struct('i',[],'j',[]);
index = struct('M',i_j,'L',i_j,'D',i_j,'U',i_j);

% Convert A to gpuArray.
A = gpuArray(A);

% Loop over the blocks to perform the factorization.
for i = 1:(s_l - 1)
    
    j = update_index_nested(s(i),j,m,n);
    
    [index_m_L,index_n_L] = GPUArray_select_local(index.L.i,index.L.j,m,n);
    [index_m_D,index_n_D] = GPUArray_select_local(index.D.i,index.D.j,m,n);
    [index_m_U,index_n_U] = GPUArray_select_local(index.U.i,index.U.j,m,n);
    
    D_inv = Calculate_D_inverse_local(A,index_m_D,index_n_D,m,n);
    
    schur_complement = (index_m_L*A*index_n_L)*D_inv*(index_m_U*A*index_n_U);
        
    A = A - schur_complement;
    
end

% Reset for the next loop.
j = 0;
A = gather(A);

% Preallocate output matrices
L_cell = cell(size(s));
D_cell = L_cell; U_cell = L_cell;

% Loop over the blocks to extract L, D, U.
for i = 1:s_l
   
    j = update_index_nested(s(i),j,m,n);
    
    L_cell{i} = fill_cells_local(A,L_cell{i},index.L.i,index.L.j);
    D_cell{i} = fill_cells_local(A,D_cell{i},index.D.i,index.D.j);
    U_cell{i} = fill_cells_local(A,U_cell{i},index.U.i,index.U.j);
    
    shift_index_nested(i,s_l);
           
end

L = collaps_cell_local(L_cell,m,n);
D = collaps_cell_local(D_cell,m,n);
U = collaps_cell_local(U_cell,m,n);

% Nested index update function
function [ j ] = update_index_nested( s_i, j, m, n )
   
    index.D.i = (1:s_i) + j;
    index.D.j = index.D.i;
    index.L.i = (index.D.i(end) + 1):m;
    index.L.j = index.D.j;
    index.U.i = index.D.i;
    index.U.j = (index.D.j(end) + 1):n;
    index.M.i = index.L.i;
    index.M.j = index.U.j;
    
    j = s_i + j;
       
end

% Nested index shift function.
function [ ] = shift_index_nested( i, s_end )

    D_cell{i}([1,2],:) =  D_cell{i}([1,2],:) + ...
                         [(index.D.i(1) - 1);(index.D.j(1) - 1)];
        
    % Only shift the index when not at the last itteration.
    if i ~= s_end
        
        L_cell{i}([1,2],:) =  L_cell{i}([1,2],:) + ...
                             [(index.L.i(1) - 1);(index.L.j(1) - 1)];
        
        U_cell{i}([1,2],:) =  U_cell{i}([1,2],:) + ...
                             [(index.U.i(1) - 1);(index.U.j(1) - 1)];
        
    end

end

end

% Local function that fills the cell array.
function [ M_cell ] = fill_cells_local( M, M_cell, index_i, index_j )
    
    % Insert COO representation of selected A entries.
    [M_cell(1,:),M_cell(2,:),M_cell(3,:)] = find(M(index_i,index_j));

end

% Local function that converts the cell array to a regular sparse matrices.
function [ M ] = collaps_cell_local( M_cell, m, n )
    
    % Strange matlab syntax that just works in this particalur case...
    M_col = [M_cell{:}];
    
    % Convert COO to CCS (MATLAB) sparse.
    M = sparse(M_col(1,:),M_col(2,:),M_col(3,:),m,n);
    
end

% Local function to facilitate slicing a sparse gpuArray.
function [ index_m, index_n ] = GPUArray_select_local( i, j, m, n )

index_m = gpuArray(sparse(i,i,1,m,m));
index_n = gpuArray(sparse(j,j,1,n,n));

end

function [ D_inv ] = Calculate_D_inverse_local( M, index_i, index_j, m, n )

[i,j,v] = find(index_i*M*index_j);

[i,j,v] = gather(i,j,v);

M_cpu = sparse((i - i(1) + 1),(j - j(1) + 1),v);

[x,y,z] = find(inv(M_cpu));

D_inv = gpuArray(sparse((x + i(1) - 1),(y + j(1) - 1),z,m,n));

end