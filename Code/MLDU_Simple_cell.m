function [ L, D, U ] = MLDU_Simple_cell( A, s )

% Creating indices structure
j = 0;
[m,n] = size(A);
s_l = length(s);
i_j = struct('i',[],'j',[]);
index = struct('M',i_j,'L',i_j,'D',i_j,'U',i_j);

% Preallocate output matrices
L_cell = cell(size(s));
D_cell = L_cell; U_cell = L_cell;

% Loop over the blocks
for i = 1:s_l
   
    j = update_index_nested(s(i),j,m,n);
    
    L_cell{i} = fill_cells_local(A,L_cell{i},index.L.i,index.L.j);
    D_cell{i} = fill_cells_local(A,D_cell{i},index.D.i,index.D.j);
    U_cell{i} = fill_cells_local(A,U_cell{i},index.U.i,index.U.j);
    
    shift_index_nested(i,s_l);
    
    schur_complement =  A(index.L.i,index.L.j)* ...
                       (A(index.D.i,index.D.j)\ ...
                        A(index.U.i,index.U.j));
    
    A(index.M.i,index.M.j) = A(index.M.i,index.M.j) - schur_complement;
           
end

L = collaps_cell_local(L_cell,m,n);
D = collaps_cell_local(D_cell,m,n);
U = collaps_cell_local(U_cell,m,n);

% Nested index update function
function [ j ] = update_index_nested( s, j, m, n )
   
    index.D.i = (1:s) + j;
    index.D.j = index.D.i;
    index.L.i = (index.D.i(end) + 1):m;
    index.L.j = index.D.j;
    index.U.i = index.D.i;
    index.U.j = (index.D.j(end) + 1):n;
    index.M.i = index.L.i;
    index.M.j = index.U.j;
    
    j = s + j;
       
end

% Nested index shift function.
function [ ] = shift_index_nested( i, s_end )

    D_cell{i}(1,:) = D_cell{i}(1,:) + (index.D.i(1) - 1);
    D_cell{i}(2,:) = D_cell{i}(2,:) + (index.D.j(1) - 1);
    
    % Only shift the index when not at the last itteration.
    if i ~= s_end
        
        L_cell{i}(1,:) = L_cell{i}(1,:) + (index.L.i(1) - 1);
        L_cell{i}(2,:) = L_cell{i}(2,:) + (index.L.j(1) - 1);
        U_cell{i}(1,:) = U_cell{i}(1,:) + (index.U.i(1) - 1);
        U_cell{i}(2,:) = U_cell{i}(2,:) + (index.U.j(1) - 1);
        
    end

end

end

% Local function that fills the cell array.
function [ A_cell ] = fill_cells_local( A, A_cell, index_i, index_j )
    
    % Insert COO representation of selected A entries.
    [A_cell(1,:),A_cell(2,:),A_cell(3,:)] = find(A(index_i,index_j));

end

% Local function that converts the cell array to a regular sparse matrices.
function [ A ] = collaps_cell_local( A_cell, m, n )
    
    % Strange matlab syntax that just works in this particalur case...
    A_col = [A_cell{:}];
    
    % Convert COO to CCS (MATLAB) sparse.
    A = sparse(A_col(1,:),A_col(2,:),A_col(3,:),m,n);
    
end
