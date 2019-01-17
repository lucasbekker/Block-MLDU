% =========================================================================
% Disclaimer:
%
% This function MLDU has been written by Lucas Bekker and presented for use
% as a prototype to Jos Maubach on 27 april 2016. Being a prototype, it is 
% not suited for production use.
% The algorithm implemented in this function has been developed by Jos
% Maubach, who also aided in certain design choises and debugging. The 
% actual implementation has been written by Lucas Bekker including certain 
% optimizations regarding runtime reduction.
%
% Questions about using this code should be directed at Jos Maubach, this
% function can not be used without explicit written permission from Jos 
% Maubach.
%
% No support shall be provided by either Lucas Bekker or Jos Maubach.
%
% Contact information:
% Lucas Bekker                      l.w.bekker@student.tue.nl
%
% =========================================================================


function [ L_MXs, D_MXs, U_MXs ] = MLDU_test5( A, block_sizes )

% Extract data from input matrix A for initialization of S.
[m,n] = size(A);
[i,j,v] = find(A);

% Create template to construct initialization of L, D, U. z serves as a way
% to preallocate memory. THIS COULD BE IMPROVED...
z = ones(m,1);
zz = zeros(m,1);
template = struct('MXs',[],'i',z,'j',z,'v',zz,'m',0,'n',0,'nz',0);

% Initialize matrices.
S = struct('MXs',[],'i',[i;z],'j',[j;z],'v',[v;zz],'m',m,'n',n,'nz',length(i));
L = struct(template); 
D = struct(template); 
U = struct(template);

% Remove unused elements from memory. This step prevents these elements
% from polluting the data in the loop, hence MATLAB is erroneous in stating 
% that these variables span multiple functions at this point. 
clear i j v m n template z zz

% Initialize the shift data used in the assembly of the matrices L, D, U. 
shift_L = [0,0];

% Memory preallocation for information about fill-in.
S_nz = zeros(length(block_sizes),1);
S_up_nz = zeros(length(block_sizes),1);

% Loop over all the blocks in the algorithm.
for k = 1:length(block_sizes)
    k
    % Calculation of Schur complement update benefits from the matrices
    % L_up, D_up and U_up to be shifted for ease of calculation.
    shift_L(2) = shift_L(1);
    shift_L(1) = shift_L(1) + block_sizes(k);
    shift_D    = [shift_L(2),shift_L(2)];
    shift_U    = [shift_L(2),shift_L(1)];
    
    % Split the matrix in it's L, D, U and S parts.
    [L_up,D_up,U_up,ind_S] = Block_Split_local(S,block_sizes(k));
    S_nz(k) = sum(ind_S);
    
    % COO representation of the Schur complement update, calculations
    % performed by MATLAB builtin inverse functionality and sparse
    % datatype.
    [i,j,v] = find(-L_up*inverse_MLDU_local(D_up)*U_up);
    S_up_nz(k) = length(v);
    m = size(L_up,1); n = size(U_up,2);
        
    % Construct the new S matrix for the next step of the algorithm and
    % join the matrices L, D and U with their 'up' counterparts. Nested
    % functions where employed to minimize the copy on write memory
    % intensive operations.
    Block_Join_S_nested(ind_S,block_sizes(k));
    Block_Join_L_nested(block_sizes(k));
    Block_Join_D_nested(block_sizes(k));
    Block_Join_U_nested(block_sizes(k));
    
end

% Extract actual values and store them in the MATLAB sparse datatype.
L_MXs = sparse(L.i,L.j,L.v,L.m,L.n);
D_MXs = sparse(D.i,D.j,D.v,D.m,D.n);
U_MXs = sparse(U.i,U.j,U.v,U.m,U.n);

% Nested function definitions.

function [ ] = Block_Join_L_nested( block_size )
    
    % Parameters.
    l  = length(L.i);
    ll = nnz(L_up);
    
    % Check if extra memory preallocation is likely required. The method
    % used to preallocate extra memory could be improved.
    if l < (ll + L.nz)
        L.i = [L.i; ones(l,1)];
        L.j = [L.j; ones(l,1)];
        L.v = [L.v; zeros(l,1)];
    end
    
    % Create index array for the location of the added values to the
    % matrix.
    index = L.nz + (1:ll);
    
    % Join the matrices by placing them consecutively in a COO
    % representation. The 'up' addition is shifted back to it's initial
    % location. Other parameters of the COO representation are also updated
    % accordingly.
    [L.i(index),L.j(index),L.v(index)] = find(L_up);
    
    L.i(index) = L.i(index) + shift_L(1);
    L.j(index) = L.j(index) + shift_L(2);
    
    L.m = L.m + block_size;
    L.n = L.n + block_size;
    L.nz = L.nz + ll;
    
end

function [ ] = Block_Join_D_nested( block_size )
    
    % Parameters.
    l  = length(D.i);
    ll = nnz(D_up);
    
    % Check if extra memory preallocation is likely required. The method
    % used to preallocate extra memory could be improved.
    if l < (ll + D.nz)
        D.i = [D.i; ones(l,1)];
        D.j = [D.j; ones(l,1)];
        D.v = [D.v; zeros(l,1)];
    end
    
    % Create index array for the location of the added values to the
    % matrix.
    index = D.nz + (1:ll);
    
    % Join the matrices by placing them consecutively in a COO
    % representation. The 'up' addition is shifted back to it's initial
    % location. Other parameters of the COO representation are also updated
    % accordingly.
    [D.i(index),D.j(index),D.v(index)] = find(D_up);
    
    D.i(index) = D.i(index) + shift_D(1);
    D.j(index) = D.j(index) + shift_D(2);
    
    D.m = D.m + block_size;
    D.n = D.n + block_size;
    D.nz = D.nz + ll;
           
end

function [ ] = Block_Join_U_nested( block_size )
    
    % Parameters.
    l  = length(U.i);
    ll = nnz(U_up);
    
    % Check if extra memory preallocation is likely required. The method
    % used to preallocate extra memory could be improved.
    if l < (ll + U.nz)
        U.i = [U.i; ones(l,1)];
        U.j = [U.j; ones(l,1)];
        U.v = [U.v; zeros(l,1)];
    end
    
    % Create index array for the location of the added values to the
    % matrix.
    index = U.nz + (1:ll);
    
    % Join the matrices by placing them consecutively in a COO
    % representation. The 'up' addition is shifted back to it's initial
    % location. Other parameters of the COO representation are also updated
    % accordingly.
    [U.i(index),U.j(index),U.v(index)] = find(U_up);
    
    U.i(index) = U.i(index) + shift_U(1);
    U.j(index) = U.j(index) + shift_U(2);
    
    U.m = U.m + block_size;
    U.n = U.n + block_size;
    U.nz = U.nz + ll;
    
end

function [ ] = Block_Join_S_nested( ind_S, block_size )
    
%     % Construction of the new S matrix for the next step in the algorithm.
%     % This takes place by joining the matrices in a similar fashion to L,
%     % D and U, but without memory preallocation. THIS COULD BE IMPROVED...
%     
%     S.i = [S.i(ind_S) - block_size; i];
%     S.j = [S.j(ind_S) - block_size; j];
%     S.v = [S.v(ind_S) ; v];
%     S.m = max(S.m - block_size,m);
%     S.n = max(S.n - block_size,n);

ind_S_clear = ~ind_S;

l = length(v);
ll = sum(ind_S_clear);
lll = l + sum(ind_S);

if l > ll
    S.i = [S.i; ones(lll,1)];
    S.j = [S.j; ones(lll,1)];
    S.v = [S.v; zeros(lll,1)];
    
    ind_S = [ind_S; false(lll,1)];
    ind_S_clear = [ind_S_clear; true(lll,1)];
    
    ll = ll + lll;
end

S.i(ind_S) = S.i(ind_S) - block_size;
S.j(ind_S) = S.j(ind_S) - block_size;

S.i(ind_S_clear) = [i; ones((ll -l),1)];
S.j(ind_S_clear) = [j; ones((ll -l),1)];
S.v(ind_S_clear) = [v; zeros((ll -l),1)];

    S.m = max(S.m - block_size,m);
    S.n = max(S.n - block_size,n);

end

end

% Local function definitions.

function [ L, D, U, ind_S ] = Block_Split_local( S, block_size )

% Construct the index vectors for the matrices L, D, U and S.
ind_L = (S.j <= block_size);
ind_U = (S.i <= block_size);
ind_D = (ind_L & ind_U);
ind_S = ~or(ind_L,ind_U);
ind_L = xor(ind_D,ind_L);
ind_U = xor(ind_D,ind_U);

% Fill L with data from S and shift where necessary.
L = sparse((S.i(ind_L) - block_size),S.j(ind_L),S.v(ind_L),...
           (S.m - block_size),block_size);
       
% Fill D with data from S.
D = sparse(S.i(ind_D),S.j(ind_D),S.v(ind_D),block_size,block_size);

% Fill L with data from S and shift where necessary.
U = sparse(S.i(ind_U),(S.j(ind_U) - block_size),S.v(ind_U),...
           block_size,(S.n - block_size));

end

function [ D_inv ] = inverse_MLDU_local( D )

% Redundant function definition at the moment, but placed here for future
% convenience where experementing with different ways to calculate the
% inverse of D might be an option.
D_inv = inv(D);

end