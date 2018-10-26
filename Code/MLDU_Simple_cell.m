function [ L, D, U ] = MLDU_Simple_cell( A, s )

% Creating indices structure
j = 0;
[m,n] = size(A);
s_l = length(s);
i_j = struct('i',[],'j',[]);
index = struct('M',i_j,'L',i_j,'D',i_j,'U',i_j);

% Preallocate output matrices
L_cell = cell(size(s));
D_cell = cell(size(s));
U_cell = cell(size(s));

% Loop over the blocks
for i = 1:s_l
   
    update_index_nested(i);
    
    fill_cells_nested(i);
    
    schur_complement =  A(index.L.i,index.L.j)* ...
                       (A(index.D.i,index.D.j)\ ...
                        A(index.U.i,index.U.j));
    
    A(index.M.i,index.M.j) = A(index.M.i,index.M.j) - schur_complement;
           
end

collaps_cells();

% Nested index function
function [ ] = update_index_nested( i )
   
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

% Nested function that fills the cell arrays of L, D and U.
function [ ] = fill_cells_nested( i )
    
    [L_cell{i}(1,:),L_cell{i}(2,:),L_cell{i}(3,:)] = ...
        find(A(index.L.i,index.L.j));
    
    [D_cell{i}(1,:),D_cell{i}(2,:),D_cell{i}(3,:)] = ...
        find(A(index.D.i,index.D.j));
    
    [U_cell{i}(1,:),U_cell{i}(2,:),U_cell{i}(3,:)] = ...
        find(A(index.U.i,index.U.j));
   
    % Shift the index
    D_cell{i}(1,:) = D_cell{i}(1,:) + (index.D.i(1) - 1);
    D_cell{i}(2,:) = D_cell{i}(2,:) + (index.D.j(1) - 1);
    
    % Only shift the index when not at the last itteration.
    if i ~= s_l
        
        L_cell{i}(1,:) = L_cell{i}(1,:) + (index.L.i(1) - 1);
        L_cell{i}(2,:) = L_cell{i}(2,:) + (index.L.j(1) - 1);
        U_cell{i}(1,:) = U_cell{i}(1,:) + (index.U.i(1) - 1);
        U_cell{i}(2,:) = U_cell{i}(2,:) + (index.U.j(1) - 1);
        
    end
    
end

% Nested function that converts the cell arrays to regular spare matrices.
function [ ] = collaps_cells( )
    
    L_col = [L_cell{:}];
    D_col = [D_cell{:}];
    U_col = [U_cell{:}];

    L = sparse(L_col(1,:),L_col(2,:),L_col(3,:),m,n);
    D = sparse(D_col(1,:),D_col(2,:),D_col(3,:));
    U = sparse(U_col(1,:),U_col(2,:),U_col(3,:),m,n);
    
end

end