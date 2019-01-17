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


function [ L_MXs, D_MXs, U_MXs ] = MLDU( A, block_sizes )

% Extract data from input matrix A for initialization of S.
[m,n] = size(A);
[i,j,v] = find(A);

% Create template to construct initialization of L, D, U. z serves as a way
% to preallocate memory. THIS COULD BE IMPROVED...
z = zeros(m,1);
template = struct('MXs',[],'i',z,'j',z,'v',z,'m',0,'n',0,'nz',0);

% Initialize matrices.
S = struct('MXs',[],'i',i,'j',j,'v',v,'m',m,'n',n);
L = struct(template); 
D = struct(template); 
U = struct(template);

% Remove unused elements from memory. This step prevents these elements
% from polluting the data in the loop, hence MATLAB is erroneous in stating 
% that these variables span multiple functions at this point. 
clear i j v m n template z

% Initialize the shift data used in the assembly of the matrices L, D, U. 
shift_L = [0,0];

% Loop over all the blocks in the algorithm.
for k = 1:length(block_sizes)
    
    % Calculation of Schur complement update benefits from the matrices
    % L_up, D_up and U_up to be shifted for ease of calculation.
    shift_L(2) = shift_L(1);
    shift_L(1) = shift_L(1) + block_sizes(k);
    shift_D    = [shift_L(2),shift_L(2)];
    shift_U    = [shift_L(2),shift_L(1)];
    
    % Split the matrix in it's L, D, U and S parts.
    [L_up,D_up,U_up,S] = Block_Split_local(S,block_sizes(k));
    
    % COO representation of the Schur complement update, calculations
    % performed by MATLAB builtin inverse functionality and sparse
    % datatype.
    [i,j,v] = find(-L_up.MXs*inverse_MLDU_local(D_up.MXs)*U_up.MXs);
    m = L_up.m; n = U_up.n;
    
    % Construct the new S matrix for the next step of the algorithm and
    % join the matrices L, D and U with their 'up' counterparts. Nested
    % functions where employed to minimize the copy on write memory
    % intensive operations.
    Block_Join_S_nested();
    Block_Join_L_nested();
    Block_Join_D_nested();
    Block_Join_U_nested();

end

% Create index arrays of actual values, the remaining values will be zero 
% and are only used as a preallocation tool.
ind_L = 1:L.nz; ind_D = 1:D.nz; ind_U = 1:U.nz;

% Extract actual values and store them in the MATLAB sparse datatype.
L_MXs = sparse(L.i(ind_L),L.j(ind_L),L.v(ind_L),L.m,L.n);
D_MXs = sparse(D.i(ind_D),D.j(ind_D),D.v(ind_D),D.m,D.n);
U_MXs = sparse(U.i(ind_U),U.j(ind_U),U.v(ind_U),U.m,U.n);

% Nested function definitions.

function [ ] = Block_Join_L_nested( )
    
    % Parameters.
    l = length(L.i);
    ll = length(L_up.i);
    
    % Check if extra memory preallocation is likely required. The method
    % used to preallocate extra memory could be improved.
    if l < (ll + L.nz)
        L.i = vertcat(L.i,zeros(l,1));
        L.j = vertcat(L.j,zeros(l,1));
        L.v = vertcat(L.v,zeros(l,1));
    end
    
    % Create index array for the location of the added values to the
    % matrix.
    index = L.nz + (1:ll);
    
    % Join the matrices by placing them consecutively in a COO
    % representation. The 'up' addition is shifted back to it's initial
    % location. Other parameters of the COO representation are also updated
    % accordingly.
    L.i(index) = (L_up.i + shift_L(1));
    L.j(index) = (L_up.j + shift_L(2));
    L.v(index) = L_up.v;
    L.m = max(L.m,L_up.m + shift_L(1));
    L.n = max(L.n,L_up.n + shift_L(2));
    L.nz = L.nz + ll;
    
end

function [ ] = Block_Join_D_nested( )
    
    % Parameters.
    l = length(D.i);
    ll = length(D_up.i);
    
    % Check if extra memory preallocation is likely required. The method
    % used to preallocate extra memory could be improved.
    if l < (ll + D.nz)
        D.i = vertcat(D.i,zeros(l,1));
        D.j = vertcat(D.j,zeros(l,1));
        D.v = vertcat(D.v,zeros(l,1));
    end
    
    % Create index array for the location of the added values to the
    % matrix.
    index = D.nz + (1:ll);
    
    % Join the matrices by placing them consecutively in a COO
    % representation. The 'up' addition is shifted back to it's initial
    % location. Other parameters of the COO representation are also updated
    % accordingly.
    D.i(index) = (D_up.i + shift_D(1));
    D.j(index) = (D_up.j + shift_D(2));
    D.v(index) = D_up.v;
    D.m = max(D.m,D_up.m + shift_D(1));
    D.n = max(D.n,D_up.n + shift_D(2));
    D.nz = D.nz + ll;
           
end

function [ ] = Block_Join_U_nested( )
    
    % Parameters.
    l = length(U.i);
    ll = length(U_up.i);
    
    % Check if extra memory preallocation is likely required. The method
    % used to preallocate extra memory could be improved.
    if l < (ll + U.nz)
        U.i = vertcat(U.i,zeros(l,1));
        U.j = vertcat(U.j,zeros(l,1));
        U.v = vertcat(U.v,zeros(l,1));
    end
    
    % Create index array for the location of the added values to the
    % matrix.
    index = U.nz + (1:ll);
    
    % Join the matrices by placing them consecutively in a COO
    % representation. The 'up' addition is shifted back to it's initial
    % location. Other parameters of the COO representation are also updated
    % accordingly.
    U.i(index) = (U_up.i + shift_U(1));
    U.j(index) = (U_up.j + shift_U(2));
    U.v(index) = U_up.v;
    U.m = max(U.m,U_up.m + shift_U(1));
    U.n = max(U.n,U_up.n + shift_U(2));
    U.nz = U.nz + ll;
    
end

function [ ] = Block_Join_S_nested( )
    
    % Construction of the new S matrix for the next step in the algorithm.
    % This takes place by joining the matrices in a similar fashion to L,
    % D and U, but without memory preallocation. THIS COULD BE IMPROVED...
    S.i = vertcat(S.i,i);
    S.j = vertcat(S.j,j);
    S.v = vertcat(S.v,v);
    S.m = max(S.m,m);
    S.n = max(S.n,n);
           
end

end

% Local function definitions.

function [ L, D, U, S ] = Block_Split_local( S, block_size )

% Construct empty structs for L, D and U.
L = struct('MXs',[],'i',[],'j',[],'v',[],'m',[],'n',[]);
D = struct(L); U = struct(L);

% Construct the index vectors for the matrices L, D, U and S.
ind_L = (S.j <= block_size);
ind_U = (S.i <= block_size);
ind_D = (ind_L & ind_U);
ind_S = ~or(ind_L,ind_U);
ind_L = xor(ind_D,ind_L);
ind_U = xor(ind_D,ind_U);

% Fill L with data from S and shift where necessary, also convert COO to
% MATLAB sparse format.
L.i = S.i(ind_L) - block_size;
L.j = S.j(ind_L);
L.v = S.v(ind_L);
L.m = S.m - block_size;
L.n = block_size;
L.MXs = sparse(L.i,L.j,L.v,L.m,L.n);

% Fill D with data from S, also convert COO to MATLAB sparse format.
D.i = S.i(ind_D);
D.j = S.j(ind_D);
D.v = S.v(ind_D);
D.m = block_size;
D.n = block_size;
D.MXs = sparse(D.i,D.j,D.v,D.m,D.n);

% Fill L with data from S and shift where necessary, also convert COO to
% MATLAB sparse format.
U.i = S.i(ind_U);
U.j = S.j(ind_U) - block_size;
U.v = S.v(ind_U);
U.m = block_size;
U.n = S.n - block_size;
U.MXs = sparse(U.i,U.j,U.v,U.m,U.n);

% Create new S.
S.i = S.i(ind_S) - block_size;
S.j = S.j(ind_S) - block_size;
S.v = S.v(ind_S);
S.m = S.m - block_size;
S.n = S.n - block_size;

end

function [ D_inv ] = inverse_MLDU_local( D )

% Redundant function definition at the moment, but placed here for future
% convenience where experementing with different ways to calculate the
% inverse of D might be an option.
D_inv = inv(D);

end