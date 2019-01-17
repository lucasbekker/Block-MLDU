function [ L_MXs, D_MXs, U_MXs ] = MLDU_test( A, block_sizes )

% Extract data from input matrix for initialization of S.
[m,n] = size(A);
[i,j,v] = find(A);

% Create template to construct initialization of L, D, U. z serves as a way
% to preallocate memory. Make more educated guess with information on
% fill-in, etc.
z(m,1) = 0;
template = struct('MXs',[],'i',z,'j',z,'v',z,'m',0,'n',0,'nz',0);

% Initialize matrices.
S = struct('MXs',[],'i',i,'j',j,'v',v,'m',m,'n',n);
L = struct(template); 
D = struct(template); 
U = struct(template);

% Remove unused elements from memory. This step prevents these elements
% from poluting the data in the loop, hence MATLAB is erroneous in stating 
% that these variables span multiple functions at this point. 
clear i j v m n template z

% Initialize the shift data used in the assembly of the matrices L, D, U. 
shift_D = [0,0];
shift_L = [block_sizes(1),0];
shift_U = [0,block_sizes(1)];

for k = 1:length(block_sizes)
    
    [L_up,D_up,U_up,S] = Block_Split_local(S,block_sizes(k));
    
    [i,j,v] = find(-L_up.MXs*inverse_MLDU_local(D_up.MXs)*U_up.MXs);
    m = L_up.m; n = U_up.n;
    
    Block_Join_S_nested();
    Block_Join_L_nested();
    Block_Join_D_nested();
    Block_Join_U_nested();

    shift_L = shift_L + block_sizes(k);
    shift_D = shift_D + block_sizes(k);
    shift_U = shift_U + block_sizes(k);
    
end

ind_L = 1:L.nz; ind_D = 1:D.nz; ind_U = 1:U.nz;

L_MXs = sparse(L.i(ind_L),L.j(ind_L),L.v(ind_L),L.m,L.n);
D_MXs = sparse(D.i(ind_D),D.j(ind_D),D.v(ind_D),D.m,D.n);
U_MXs = sparse(U.i(ind_U),U.j(ind_U),U.v(ind_U),U.m,U.n);

function [ ] = Block_Join_L_nested( )
    
    l = length(L.i);
    ll = length(L_up.i);
    
    if l < (ll + L.nz)
        L.i(2*l) = 0;
        L.j(2*l) = 0;
        L.v(2*l) = 0;
    end
    
    index = L.nz + (1:ll);
    
    L.i(index) = (L_up.i + shift_L(1));
    L.j(index) = (L_up.j + shift_L(2));
    L.v(index) = L_up.v;
    L.m = max(L.m,L_up.m + shift_L(1));
    L.n = max(L.n,L_up.n + shift_L(2));
    L.nz = L.nz + ll;
            
end
function [ ] = Block_Join_D_nested( )
    
    l = length(D.i);
    ll = length(D_up.i);
    
    if l < (ll + D.nz)
        D.i(2*l) = 0;
        D.j(2*l) = 0;
        D.v(2*l) = 0;
    end
    
    index = D.nz + (1:ll);
    
    D.i(index) = (D_up.i + shift_D(1));
    D.j(index) = (D_up.j + shift_D(2));
    D.v(index) = D_up.v;
    D.m = max(D.m,D_up.m + shift_D(1));
    D.n = max(D.n,D_up.n + shift_D(2));
    D.nz = D.nz + ll;
           
end
function [ ] = Block_Join_U_nested( )
        
    l = length(U.i);
    ll = length(U_up.i);
    
    if l < (ll + U.nz)
        U.i(2*l) = 0;
        U.j(2*l) = 0;
        U.v(2*l) = 0;
    end
    
    index = U.nz + (1:ll);
    
    U.i(index) = (U_up.i + shift_U(1));
    U.j(index) = (U_up.j + shift_U(2));
    U.v(index) = U_up.v;
    U.m = max(U.m,U_up.m + shift_U(1));
    U.n = max(U.n,U_up.n + shift_U(2));
    U.nz = U.nz + ll;
    
end

function [ ] = Block_Join_S_nested( )
    
    S.i = vertcat(S.i,i);
    S.j = vertcat(S.j,j);
    S.v = vertcat(S.v,v);
    S.m = max(S.m,m);
    S.n = max(S.n,n);
    
end

end

function [ L, D, U, S ] = Block_Split_local( S, block_size )

L = struct('MXs',[],'i',[],'j',[],'v',[],'m',[],'n',[]);
D = struct(L); U = struct(L);

ind_L = logical(bsxfun(@le,S.j,block_size));
ind_U = logical(bsxfun(@le,S.i,block_size));
ind_D = bitand(ind_L,ind_U);
ind_S = ~bitor(ind_L,ind_U);
ind_L = bitxor(ind_D,ind_L);
ind_U = bitxor(ind_D,ind_U);

L.i = S.i(ind_L) - block_size;
L.j = S.j(ind_L);
L.v = S.v(ind_L);
L.m = S.m - block_size;
L.n = block_size;
L.MXs = sparse(L.i,L.j,L.v,L.m,L.n);

D.i = S.i(ind_D);
D.j = S.j(ind_D);
D.v = S.v(ind_D);
D.m = block_size;
D.n = block_size;
D.MXs = sparse(D.i,D.j,D.v,D.m,D.n);

U.i = S.i(ind_U);
U.j = S.j(ind_U) - block_size;
U.v = S.v(ind_U);
U.m = block_size;
U.n = S.n - block_size;
U.MXs = sparse(U.i,U.j,U.v,U.m,U.n);

S.i = S.i(ind_S) - block_size;
S.j = S.j(ind_S) - block_size;
S.v = S.v(ind_S);
S.m = S.m - block_size;
S.n = S.n - block_size;

end

function [ D_inv ] = inverse_MLDU_local( D )

D_inv = D^(-1);

end